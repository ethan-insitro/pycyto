"""Per-probe and run-level QC metrics for a ``cyto workflow gex`` run."""

import logging
import os
from typing import Any

import anndata as ad
import numpy as np
import polars as pl
import scipy.sparse as sp

from .parse import UNMAPPED_LABELS, load_json, read_barcode_stats

logger = logging.getLogger("pycyto.qc")

# log10 bin edges (width 0.05) for the histograms embedded in the report
LOG_BINS = np.round(np.arange(0, 6.05, 0.05), 2)


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


def rank_curve(umis_desc: np.ndarray, is_cell_desc: np.ndarray, n_points: int = 300) -> list:
    """Barcode-rank curve downsampled to ~n_points log-spaced ranks.

    Each point is ``[rank, umis, fraction of barcodes in the segment that are cells]``.
    """
    n = len(umis_desc)
    if n == 0:
        return []
    ranks = np.unique(np.clip(np.round(np.logspace(0, np.log10(n), n_points)), 1, n)).astype(int)
    prev = np.concatenate([[0], ranks[:-1]])
    csum = np.concatenate([[0], np.cumsum(is_cell_desc)])
    frac = (csum[ranks] - csum[prev]) / (ranks - prev)
    return [[int(r), int(u), round(float(f), 3)] for r, u, f in zip(ranks, umis_desc[ranks - 1], frac)]


def log_hist(values: np.ndarray) -> list[int] | None:
    if len(values) == 0:
        return None
    counts, _ = np.histogram(np.log10(np.maximum(values, 1)), bins=np.append(LOG_BINS, 6.05))
    return counts.tolist()


def unmapped_reasons(mapping: dict) -> list[dict[str, Any]]:
    """Unmapped-read categories, largest first. A read can fail more than one check."""
    unmapped = mapping["unmapped"]
    rows = [
        {
            "reason": key,
            "label": UNMAPPED_LABELS.get(key, key.replace("_", " ")),
            "reads": reads,
            "frac_of_reads": _div(reads, mapping["total_reads"]),
            "frac_of_unmapped": unmapped[f"{key}_frac"],
        }
        for key, reads in unmapped.items()
        if not key.endswith("_frac")
    ]
    return sorted(rows, key=lambda r: -r["reads"])


def process_probe(cyto_outdir: str, probe: str) -> dict[str, Any]:
    """Metrics, plot data and per-cell arrays for one probe barcode.

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
    df = df.with_columns(is_cell=pl.col("n_genes").is_not_null()).sort("n_umis", descending=True)
    in_cells = df.filter("is_cell")
    cell_umis = in_cells["n_umis"].to_numpy()
    cell_genes = in_cells["n_genes"].to_numpy()
    mapped = df["n_reads"].sum()
    return {
        "rec": {
            "probe": probe,
            "mapped_reads": mapped,
            "umis": df["n_umis"].sum(),
            "cells": in_cells.height,
            "reads_in_cells": in_cells["n_reads"].sum(),
            "frac_reads_in_cells": _div(in_cells["n_reads"].sum(), mapped),
            "median_umis_per_cell": in_cells["n_umis"].median(),
        },
        "plots": {
            "curve": rank_curve(df["n_umis"].to_numpy(), df["is_cell"].to_numpy()),
            "umi_hist": log_hist(cell_umis),
            "gene_hist": log_hist(cell_genes),
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
    reasons = unmapped_reasons(mapping)

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
        "top_unmapped_reason": reasons[0]["label"] if reasons else None,
        "failed_umi_qual_of_total": _div(mapping["unmapped"]["failed_umi_qual"], mapping["total_reads"]),
        "probe_barcodes_in_library": lib["probe"]["total_elem"],
        "probe_barcodes_with_reads": len(recs),
        "seq_saturation": 1 - umis / mapped if mapped else None,
        "umi_corrected_frac": _div(corrected, total_umis),
        # cells (cyto's filtered h5ad)
        "estimated_cells": n_cells,
        "probe_barcodes_with_cells": called.height,
        "n_probes_without_cells": probes.height - called.height,
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


def pooled_plots(results: list[dict]) -> dict[str, Any]:
    """Run-wide histograms: UMIs and genes per cell across all probe barcodes."""
    return {
        "umi_hist": log_hist(np.concatenate([r["cell_umis"] for r in results])),
        "gene_hist": log_hist(np.concatenate([r["cell_genes"] for r in results])),
    }
