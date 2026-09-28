"""Cell Ranger-style QC report for a single ``cyto workflow gex`` output directory.

Entry point: :func:`build_report` (CLI: ``pycyto qc``).
"""

import datetime as dt
import logging
import multiprocessing as mp
import os
from functools import partial

import numpy as np

from .alerts import build_alerts
from .metrics import LOG_BINS, log_hist, process_probe, summarize
from .parse import discover_probes, load_json, parse_cyto_log, parse_timings
from .render import render_html, write_csvs

__all__ = ["build_report", "collect"]

logger = logging.getLogger("pycyto.qc")

# per-probe columns that are only needed internally
_INTERNAL_COLS = ("n_vars", "umi_corrected", "umi_total")


def _pycyto_version() -> str:
    try:
        from importlib.metadata import version

        return version("pycyto")
    except Exception:
        return "unknown"


def _init_worker(verbose: bool = False) -> None:
    """Initialize logging in each spawned worker (mirrors aggregate.init_worker)."""
    import sys

    worker_logger = logging.getLogger("pycyto.qc")
    worker_logger.handlers.clear()
    handler = logging.StreamHandler(sys.stderr)
    handler.setFormatter(
        logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    )
    worker_logger.addHandler(handler)
    worker_logger.setLevel(logging.DEBUG if verbose else logging.INFO)
    worker_logger.propagate = False


def collect(
    cyto_outdir: str, threads: int = -1, verbose: bool = False, title: str | None = None
) -> dict:
    """Compute every metric and plot input for the report; returns the report payload."""
    if not os.path.isdir(os.path.join(cyto_outdir, "stats")):
        raise FileNotFoundError(f"{cyto_outdir} doesn't look like a cyto output directory (no stats/)")
    probes = discover_probes(cyto_outdir)
    if not probes:
        raise FileNotFoundError(f"No stats/reads/*.reads.tsv.zst found under {cyto_outdir}")

    stats = os.path.join(cyto_outdir, "stats")
    mapping = load_json(os.path.join(stats, "mapping_map.json"), {}) or {}
    lib = load_json(os.path.join(stats, "mapping_lib.json"), []) or []
    run = load_json(os.path.join(stats, "mapping_run.json"), []) or []
    loginfo = parse_cyto_log(os.path.join(cyto_outdir, "cyto.log"))
    timings = parse_timings(os.path.join(cyto_outdir, ".timings.tsv"))

    if threads == -1:
        threads = mp.cpu_count()
    threads = max(1, min(threads, len(probes)))
    logger.info(f"Computing QC for {len(probes)} probe barcodes using {threads} threads")

    func = partial(process_probe, cyto_outdir)
    if threads == 1:
        results = [func(p) for p in probes]
    else:
        ctx = mp.get_context("spawn")
        with ctx.Pool(threads, initializer=partial(_init_worker, verbose=verbose)) as pool:
            results = pool.map(func, probes)

    recs = [r["rec"] for r in results]
    extras = [r["extras"] for r in results]
    total_mapped = sum(r["mapped_reads"] for r in recs)
    for r in recs:
        r["frac_of_mapped_reads"] = r["mapped_reads"] / total_mapped if total_mapped else None

    cell_umis = np.concatenate([e["cell_umis"] for e in extras])
    cell_genes = np.concatenate([e["cell_genes"] for e in extras])
    detected = [e["detected"] for e in extras if e["detected"] is not None]

    summary = summarize(
        recs, cell_umis, cell_genes, detected, mapping, lib, loginfo, timings, len(run), cyto_outdir
    )
    public = [{k: v for k, v in r.items() if k not in _INTERNAL_COLS} for r in recs]
    return {
        "title": title or os.path.basename(os.path.abspath(cyto_outdir)),
        "generated": dt.datetime.now().isoformat(sep=" ", timespec="seconds"),
        "version": _pycyto_version(),
        "summary": summary,
        "alerts": build_alerts(summary, recs),
        "unmapped": mapping.get("unmapped", {}) or {},
        "library": lib,
        "run": run,
        "timings": timings,
        "probes": public,
        "extras": {
            r["probe"]: {"curve": e["curve"], "umi_hist": e["umi_hist"], "gene_hist": e["gene_hist"]}
            for r, e in zip(recs, extras)
        },
        "pooled": {
            "umi_hist": _hist_or_none(cell_umis),
            "gene_hist": _hist_or_none(cell_genes),
        },
        "log_bins": LOG_BINS.tolist(),
    }


def _hist_or_none(values: np.ndarray) -> list[int] | None:
    return log_hist(values) if len(values) else None


def build_report(
    cyto_outdir: str,
    output: str | None = None,
    title: str | None = None,
    write_csv: bool = True,
    threads: int = -1,
    verbose: bool = False,
) -> str:
    """Write the HTML report (and optionally metric CSVs). Returns the HTML path."""
    payload = collect(cyto_outdir, threads=threads, verbose=verbose, title=title)
    output = output or os.path.join(cyto_outdir, "qc_report.html")
    with open(output, "w", encoding="utf-8") as fh:
        fh.write(render_html(payload))
    logger.info(f"Wrote QC report: {output}")
    if write_csv:
        stem = os.path.splitext(output)[0]
        paths = write_csvs(payload["summary"], payload["probes"], stem)
        logger.info(f"Wrote metrics: {', '.join(paths)}")
    for a in payload["alerts"]:
        if a["level"] != "ok":
            logger.warning(f"{a['title']}: {a['detail']}")
    return output
