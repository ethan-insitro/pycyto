"""Per-probe and run-level QC metrics for a ``cyto workflow gex`` run."""

import logging
import os
from typing import Any

import anndata as ad
import numpy as np
import polars as pl
import scipy.sparse as sp

from .parse import load_json, read_barcode_stats

logger = logging.getLogger("pycyto.qc")


def _div(a, b) -> float | None:
    return a / b if a is not None and b else None


def read_counts(path: str, chunk_rows: int = 10_000) -> tuple[pl.DataFrame, np.ndarray, list[str]]:
    """Scan a cyto count h5ad.

    Returns a ``barcode, n_features`` frame (features with nonzero counts per barcode),
    the per-feature count totals, and the feature names. ``X`` is read in row chunks
    from a backed AnnData so memory stays bounded.
    """
    adata = ad.read_h5ad(path, backed="r")
    try:
        n_features = np.zeros(adata.n_obs, dtype=np.int64)
        totals = np.zeros(adata.n_vars, dtype=np.float64)
        for start in range(0, adata.n_obs, chunk_rows):
            block = sp.csr_matrix(adata.X[start : start + chunk_rows])
            block.eliminate_zeros()
            n_features[start : start + block.shape[0]] = np.diff(block.indptr)
            totals += np.asarray(block.sum(axis=0)).ravel()
        # cyto names cells ``<barcode>-<probe>`` (e.g. ``ACGT...-A-A02``)
        barcodes = [name.split("-", 1)[0] for name in adata.obs_names]
        names = adata.var_names.tolist()
    finally:
        adata.file.close()
    return pl.DataFrame({"barcode": barcodes, "n_features": n_features}), totals, names


def process_probe(cyto_outdir: str, probe: str) -> dict[str, Any]:
    """Metrics and per-cell arrays for one probe barcode.

    Cells are exactly the barcodes in cyto's ``counts/<probe>.filt.h5ad``; probe barcodes
    without that file have no cells. ``rec`` is the row shown in the report's probe table;
    ``umi_counts`` is the ``(corrected, total)`` pair from ``stats/umi/<probe>.umi.json``,
    only used by :func:`summarize`.
    """
    stats = os.path.join(cyto_outdir, "stats")
    df = read_barcode_stats(os.path.join(stats, "reads", f"{probe}.reads.tsv.zst"))
    umi = load_json(os.path.join(stats, "umi", f"{probe}.umi.json"))

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
    df = df.with_columns(is_cell=pl.col("n_genes").is_not_null())
    in_cells = df.filter("is_cell")
    cell_umis = in_cells["n_umis"].to_numpy()
    cell_genes = in_cells["n_genes"].to_numpy()
    return {
        "rec": {
            "probe": probe,
            "mapped_reads": df["n_reads"].sum(),
            "umis": df["n_umis"].sum(),
            "cells": in_cells.height,
            "reads_in_cells": in_cells["n_reads"].sum(),
        },
        "cell_umis": cell_umis,
        "cell_genes": cell_genes,
        "detected": detected,
        "umi_counts": (umi["corrected"], umi["total"]),
    }


def summarize(results: list[dict], meta: dict, cyto_outdir: str) -> dict[str, Any]:
    """Run-level metrics (the numbers behind the report's Summary tab)."""
    recs = [r["rec"] for r in results]
    mapped = sum(r["mapped_reads"] for r in recs)
    umis = sum(r["umis"] for r in recs)
    corrected = sum(r["umi_counts"][0] for r in results)
    total_umis = sum(r["umi_counts"][1] for r in results)
    mapping = meta["mapping"]
    lib = {d["name"]: d for d in meta["library"]}

    probes = pl.DataFrame(recs, infer_schema_length=None)
    called = probes.filter(pl.col("cells") > 0)
    n_cells = called["cells"].sum()
    cell_umis = np.concatenate([r["cell_umis"] for r in results])
    cell_genes = np.concatenate([r["cell_genes"] for r in results])
    detected = [r["detected"] for r in results if r["detected"] is not None]
    return {
        "cyto_outdir": os.path.abspath(cyto_outdir),
        "total_reads": mapping["total_reads"],
        "mapped_reads": mapping["mapped_reads"],
        "mapped_reads_frac": mapping["mapped_reads_frac"],
        "probe_barcodes_in_library": lib["probe"]["total_elem"],
        "probe_barcodes_with_reads": len(recs),
        "seq_saturation": 1 - umis / mapped if mapped else None,
        "umi_corrected_frac": _div(corrected, total_umis),
        # cells (cyto's filtered h5ad)
        "estimated_cells": n_cells,
        "probe_barcodes_with_cells": called.height,
        "mean_reads_per_cell": _div(mapping["total_reads"], n_cells),
        "mean_mapped_reads_per_cell": _div(mapping["mapped_reads"], n_cells),
        "median_umis_per_cell": float(np.median(cell_umis)) if len(cell_umis) else None,
        "median_genes_per_cell": float(np.median(cell_genes)) if len(cell_genes) else None,
        "total_genes_detected": int(np.logical_or.reduce(detected).sum()) if detected else None,
        "genes_in_reference": len(detected[0]) if detected else lib["gex"]["total_aggr"],
        "frac_reads_in_cells": _div(probes["reads_in_cells"].sum(), mapped),
        "background_probe_read_frac": _div(probes.filter(pl.col("cells") == 0)["mapped_reads"].sum(), mapped),
        "cells_median_per_probe": called["cells"].median(),
    }
