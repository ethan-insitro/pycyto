"""GEX (``cyto workflow gex``) metrics.

Cells are exactly the barcodes in cyto's ``counts/<probe>.filt.h5ad``; probe barcodes
without that file have no cells.
"""

import logging
import os
from typing import Any

import numpy as np
import polars as pl

from .metrics import _div, log_hist, probe_basics, rank_curve, read_counts, run_summary

logger = logging.getLogger("pycyto.qc")

FEATURE = "gene probe"  # what an unmapped "missing_feature" read failed to match


def process_probe(cyto_outdir: str, probe: str) -> dict[str, Any]:
    """Metrics, plot data and per-cell arrays for one probe barcode."""
    df, rec, umi_counts = probe_basics(cyto_outdir, probe)

    filt = os.path.join(cyto_outdir, "counts", f"{probe}.filt.h5ad")
    detected = None
    if os.path.exists(filt):
        cells, totals, _ = read_counts(filt)
        detected = totals > 0
        df = df.join(cells.rename({"n_features": "n_genes"}), on="barcode", how="left")
        if (missing := cells.height - df["n_genes"].count()) > 0:
            logger.warning(f"[{probe}] - {missing} filtered barcodes missing from reads stats")
    else:
        df = df.with_columns(n_genes=pl.lit(None, dtype=pl.Int64))
    df = df.with_columns(is_cell=pl.col("n_genes").is_not_null()).sort("n_umis", descending=True)
    in_cells = df.filter("is_cell")

    n_cells, mapped = in_cells.height, rec["mapped_reads"]
    rec |= {
        "cells": n_cells,
        "has_filtered_h5ad": detected is not None,
        "reads_in_cells": in_cells["n_reads"].sum(),
        "frac_reads_in_cells": _div(in_cells["n_reads"].sum(), mapped),
        "frac_umis_in_cells": _div(in_cells["n_umis"].sum(), rec["umis"]),
        "mean_reads_per_cell": _div(mapped, n_cells),
        "median_umis_per_cell": in_cells["n_umis"].median(),
        "median_genes_per_cell": in_cells["n_genes"].median(),
        "total_genes_detected": int(detected.sum()) if detected is not None else None,
    }
    cell_umis = in_cells["n_umis"].to_numpy()
    cell_genes = in_cells["n_genes"].to_numpy()
    return {
        "rec": rec,
        "plots": {
            "curve": rank_curve(df["n_umis"].to_numpy(), df["is_cell"].to_numpy()),
            "umi_hist": log_hist(cell_umis),
            "gene_hist": log_hist(cell_genes),
        },
        "cell_umis": cell_umis,
        "cell_genes": cell_genes,
        "detected": detected,
        "umi_counts": umi_counts,
    }


def summarize(results: list[dict], meta: dict, cyto_outdir: str) -> dict[str, Any]:
    """Run-level metrics (the numbers behind the report's Summary tab)."""
    probes = pl.DataFrame([r["rec"] for r in results], infer_schema_length=None)
    called = probes.filter(pl.col("cells") > 0)
    lib = {d["name"]: d for d in meta["library"]}
    cell_umis = np.concatenate([r["cell_umis"] for r in results])
    cell_genes = np.concatenate([r["cell_genes"] for r in results])
    detected = [r["detected"] for r in results if r["detected"] is not None]
    total_reads, mapped_reads = meta["mapping"]["total_reads"], meta["mapping"]["mapped_reads"]
    n_cells, probe_mapped = called["cells"].sum(), probes["mapped_reads"].sum()
    return run_summary(results, meta, cyto_outdir, FEATURE) | {
        "estimated_cells": n_cells,
        "probe_barcodes_with_cells": called.height,
        "n_probes_without_cells": probes.height - called.height,
        "mean_reads_per_cell": _div(total_reads, n_cells),
        "mean_mapped_reads_per_cell": _div(mapped_reads, n_cells),
        "median_umis_per_cell": float(np.median(cell_umis)) if len(cell_umis) else None,
        "median_genes_per_cell": float(np.median(cell_genes)) if len(cell_genes) else None,
        "total_genes_detected": int(np.logical_or.reduce(detected).sum()) if detected else None,
        "genes_in_reference": len(detected[0]) if detected else lib["gex"]["total_aggr"],
        "frac_reads_in_cells": _div(probes["reads_in_cells"].sum(), probe_mapped),
        "background_probe_read_frac": _div(
            probes.filter(pl.col("cells") == 0)["mapped_reads"].sum(), probe_mapped
        ),
        "cells_median_per_probe": called["cells"].median(),
    }


def pooled_plots(results: list[dict]) -> dict[str, Any]:
    """Run-wide histograms: UMIs and genes per cell across all probe barcodes."""
    return {
        "umi_hist": log_hist(np.concatenate([r["cell_umis"] for r in results])),
        "gene_hist": log_hist(np.concatenate([r["cell_genes"] for r in results])),
    }
