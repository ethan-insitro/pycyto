"""Cell Ranger-style QC reports for a single cyto output directory (``cyto workflow gex``).

Entry point: :func:`build_report` (CLI: ``pycyto qc``).
"""

import datetime as dt
import logging
import os
from importlib.metadata import version

from .metrics import process_probe, summarize
from .parse import detect_workflow, discover_probes, load_run_metadata
from .render import render_html, write_csvs

__all__ = ["build_report", "collect"]

logger = logging.getLogger("pycyto.qc")


def collect(cyto_outdir: str, title: str | None = None) -> dict:
    """Compute every metric for the report; returns the report payload."""
    probes = discover_probes(cyto_outdir)
    meta = load_run_metadata(cyto_outdir)
    workflow = detect_workflow(meta)
    logger.info(f"Computing QC for a cyto {workflow} run with {len(probes)} probe barcodes")
    results = [process_probe(cyto_outdir, probe) for probe in probes]
    return {
        "workflow": workflow,
        "title": title or os.path.basename(os.path.abspath(cyto_outdir)),
        "generated": dt.datetime.now().isoformat(sep=" ", timespec="seconds"),
        "version": version("pycyto"),
        "summary": summarize(results, meta, cyto_outdir),
        "probes": [r["rec"] for r in results],
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
