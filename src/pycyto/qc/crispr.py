"""CRISPR (``cyto workflow crispr``) metrics.

A CRISPR run has no cell calls of its own: ``counts/<probe>.h5ad`` holds guide UMIs for
every barcode. The report therefore describes guide capture per probe barcode and how
evenly the guide library is covered across the run. Guide assignment
(``assignments/``) is not reported yet.
"""

import os
from typing import Any

import numpy as np

from .metrics import _div, log_hist, probe_basics, rank_curve, read_counts, run_summary

FEATURE = "guide"  # what an unmapped "missing_feature" read failed to match
TOP_GUIDES = 10  # most abundant guides listed in the report


def process_probe(cyto_outdir: str, probe: str) -> dict[str, Any]:
    """Metrics, plot data and per-guide UMI totals for one probe barcode."""
    df, rec, umi_counts = probe_basics(cyto_outdir, probe)
    _, guide_umis, guides = read_counts(os.path.join(cyto_outdir, "counts", f"{probe}.h5ad"))
    rec["guides_detected"] = int((guide_umis > 0).sum())
    umis_desc = df.sort("n_umis", descending=True)["n_umis"].to_numpy()
    return {
        "rec": rec,
        "plots": {
            "curve": rank_curve(umis_desc, np.zeros(len(umis_desc), dtype=bool)),
            "guide_hist": log_hist(guide_umis[guide_umis > 0]),
        },
        "guide_umis": guide_umis,
        "guides": guides,
        "umi_counts": umi_counts,
    }


def _guide_totals(results: list[dict]) -> tuple[np.ndarray, list[str]]:
    """UMIs per guide summed over probe barcodes (all count files share one guide list)."""
    return np.sum([r["guide_umis"] for r in results], axis=0), results[0]["guides"]


def summarize(results: list[dict], meta: dict, cyto_outdir: str) -> dict[str, Any]:
    """Run-level metrics (the numbers behind the report's Summary tab)."""
    totals, _ = _guide_totals(results)
    p10, p90 = np.percentile(totals, [10, 90])
    return run_summary(results, meta, cyto_outdir, FEATURE) | {
        "guides_in_library": len(totals),
        "guides_detected": int((totals > 0).sum()),
        "frac_guides_detected": float((totals > 0).mean()),
        "guide_umis": int(totals.sum()),
        "median_umis_per_guide": float(np.median(totals)),
        # 90th / 10th percentile of UMIs per guide: a standard library-evenness measure
        "guide_skew_ratio": _div(float(p90), float(p10)) if p10 else None,
        "mean_reads_per_probe": meta["mapping"]["mapped_reads"] / len(results),
    }


def pooled_plots(results: list[dict]) -> dict[str, Any]:
    """Run-wide guide coverage: UMIs-per-guide histogram and the most abundant guides."""
    totals, guides = _guide_totals(results)
    order = np.argsort(-totals, kind="stable")[:TOP_GUIDES]
    total = totals.sum()
    return {
        "guide_hist": log_hist(totals[totals > 0]),
        "top_guides": [
            {"guide": guides[i], "umis": int(totals[i]), "frac": _div(float(totals[i]), float(total))}
            for i in order
            if totals[i] > 0
        ],
    }
