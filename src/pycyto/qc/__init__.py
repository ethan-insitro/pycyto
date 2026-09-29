"""Cell Ranger-style QC reports for a single cyto output directory (``cyto workflow gex``).

Entry point: :func:`build_report` (CLI: ``pycyto qc``). The workflow is detected from the
output directory; each workflow module (:mod:`.gex`) provides ``process_probe``
and ``summarize``, and everything else is shared.
"""

import datetime as dt
import logging
import os
from importlib.metadata import version

import polars as pl

from . import gex
from .parse import detect_workflow, discover_probes, load_run_metadata
from .render import render_html, write_csvs

__all__ = ["build_report", "collect"]

logger = logging.getLogger("pycyto.qc")

WORKFLOWS = {"gex": gex}


def collect(cyto_outdir: str, title: str | None = None) -> dict:
    """Compute every metric for the report; returns the report payload."""
    probes = discover_probes(cyto_outdir)
    meta = load_run_metadata(cyto_outdir)
    workflow = detect_workflow(meta)
    wf = WORKFLOWS[workflow]
    logger.info(f"Computing QC for a cyto {workflow} run with {len(probes)} probe barcodes")
    results = [wf.process_probe(cyto_outdir, probe) for probe in probes]

    table = pl.DataFrame([r["rec"] for r in results], infer_schema_length=None)
    return {
        "workflow": workflow,
        "title": title or os.path.basename(os.path.abspath(cyto_outdir)),
        "generated": dt.datetime.now().isoformat(sep=" ", timespec="seconds"),
        "version": version("pycyto"),
        "summary": wf.summarize(results, meta, cyto_outdir),
        "probes": table.to_dicts(),
    }


def build_report(
    cyto_outdir: str,
    output: str | None = None,
    title: str | None = None,
    write_csv: bool = True,
) -> str:
    """Write the HTML report (and optionally the metric CSVs). Returns the HTML path."""
    payload = collect(cyto_outdir, title=title)
    output = output or os.path.join(cyto_outdir, "qc_report.html")
    with open(output, "w", encoding="utf-8") as fh:
        fh.write(render_html(payload))
    logger.info(f"Wrote QC report: {output}")
    if write_csv:
        paths = write_csvs(payload, os.path.splitext(output)[0])
        logger.info(f"Wrote metrics: {', '.join(paths)}")
    return output
