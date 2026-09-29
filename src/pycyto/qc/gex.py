"""GEX (``cyto workflow gex``) metrics."""

from typing import Any

from .metrics import probe_basics, run_summary


def process_probe(cyto_outdir: str, probe: str) -> dict[str, Any]:
    """Metrics for one probe barcode."""
    _, rec, umi_counts = probe_basics(cyto_outdir, probe)
    return {"rec": rec, "umi_counts": umi_counts}


def summarize(results: list[dict], meta: dict, cyto_outdir: str) -> dict[str, Any]:
    """Run-level metrics (the numbers behind the report's Summary tab)."""
    return run_summary(results, meta, cyto_outdir)
