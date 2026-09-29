"""Cell Ranger-style QC reports for a single cyto output directory (``cyto workflow gex``).

Entry point: :func:`build_report` (CLI: ``pycyto qc``). The workflow is detected from the
output directory; each workflow module (:mod:`.gex`) provides ``FEATURE``, ``process_probe``,
``summarize`` and ``pooled_plots``, and everything else is shared.
"""

import datetime as dt
import logging
import multiprocessing as mp
import os
import sys
from functools import partial
from importlib.metadata import version

import polars as pl

from . import gex
from .alerts import build_alerts
from .metrics import LOG_BINS
from .parse import detect_workflow, discover_probes, load_run_metadata
from .render import write_csvs

__all__ = ["build_report", "collect"]

logger = logging.getLogger("pycyto.qc")

WORKFLOWS = {"gex": gex}


def _init_worker(verbose: bool = False) -> None:
    """Log from spawned workers to stderr (their ``pycyto`` logger is unconfigured)."""
    logging.basicConfig(
        stream=sys.stderr,
        level=logging.DEBUG if verbose else logging.INFO,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    )


def _run(process_probe, cyto_outdir: str, probes: list[str], threads: int, verbose: bool) -> list[dict]:
    threads = max(1, min(mp.cpu_count() if threads == -1 else threads, len(probes)))
    logger.info(f"Computing QC for {len(probes)} probe barcodes using {threads} threads")
    func = partial(process_probe, cyto_outdir)
    if threads == 1:
        return [func(p) for p in probes]
    ctx = mp.get_context("spawn")
    with ctx.Pool(threads, initializer=_init_worker, initargs=(verbose,)) as pool:
        return pool.map(func, probes)


def collect(
    cyto_outdir: str, threads: int = -1, verbose: bool = False, title: str | None = None
) -> dict:
    """Compute every metric and plot input for the report; returns the report payload."""
    probes = discover_probes(cyto_outdir)
    if not probes:
        raise FileNotFoundError(
            f"{cyto_outdir} doesn't look like a cyto output directory "
            "(no stats/reads/*.reads.tsv.zst)"
        )
    meta = load_run_metadata(cyto_outdir)
    workflow = detect_workflow(meta)
    wf = WORKFLOWS[workflow]
    logger.info(f"Detected a cyto {workflow} run")
    results = _run(wf.process_probe, cyto_outdir, probes, threads, verbose)

    table = pl.DataFrame([r["rec"] for r in results], infer_schema_length=None)

    summary = wf.summarize(results, meta, cyto_outdir)
    return {
        "workflow": workflow,
        "title": title or os.path.basename(os.path.abspath(cyto_outdir)),
        "generated": dt.datetime.now().isoformat(sep=" ", timespec="seconds"),
        "version": version("pycyto"),
        "summary": summary,
        "alerts": build_alerts(workflow, summary, table),
        "probes": table.to_dicts(),
        "plots": {r["rec"]["probe"]: r["plots"] for r in results},
        "pooled": wf.pooled_plots(results),
        "log_bins": LOG_BINS.tolist(),
    }


def build_report(
    cyto_outdir: str, output: str | None = None, threads: int = -1, verbose: bool = False
) -> str:
    """Write the run-level metrics to ``<output>_metrics_summary.csv`` and log any alerts.

    ``output`` is a path stem (default ``<cyto_outdir>/qc_report``); returns it.
    """
    payload = collect(cyto_outdir, threads=threads, verbose=verbose)
    stem = output or os.path.join(cyto_outdir, "qc_report")
    logger.info(f"Wrote metrics: {', '.join(write_csvs(payload, stem))}")
    for alert in payload["alerts"]:
        if alert["level"] != "ok":
            logger.warning(f"{alert['title']}: {alert['detail']}")
    return stem
