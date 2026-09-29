"""QC building blocks shared by every cyto workflow (GEX, CRISPR).

Workflow-specific metrics live in :mod:`pycyto.qc.gex` and :mod:`pycyto.qc.crispr`.
"""

import os
from typing import Any

import anndata as ad
import numpy as np
import polars as pl
import scipy.sparse as sp

from .parse import UNMAPPED_LABELS, load_json, read_barcode_stats

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


def probe_basics(cyto_outdir: str, probe: str) -> tuple[pl.DataFrame, dict[str, Any], tuple[int, int]]:
    """Per-barcode reads stats plus the metrics every workflow reports for a probe barcode.

    Returns the ``barcode, n_umis, n_reads`` frame, the shared metrics and the
    ``(corrected, total)`` UMI counts from ``stats/umi/<probe>.umi.json``.
    """
    stats = os.path.join(cyto_outdir, "stats")
    df = read_barcode_stats(os.path.join(stats, "reads", f"{probe}.reads.tsv.zst"))
    umi_stats = load_json(os.path.join(stats, "umi", f"{probe}.umi.json"))
    mapped, umis = df["n_reads"].sum(), df["n_umis"].sum()
    rec = {
        "probe": probe,
        "n_barcodes": df.height,
        "mapped_reads": mapped,
        "umis": umis,
        "seq_saturation": 1 - umis / mapped if mapped else None,
        "umi_corrected_frac": umi_stats["fraction_corrected"],
    }
    return df, rec, (umi_stats["corrected"], umi_stats["total"])


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


def unmapped_reasons(mapping: dict, feature: str = "gene probe") -> list[dict[str, Any]]:
    """Unmapped-read categories, largest first. A read can fail more than one check.

    ``feature`` names what ``missing_feature`` means for the workflow (gene probe, guide).
    """
    unmapped = mapping["unmapped"]
    labels = UNMAPPED_LABELS | {"missing_feature": f"No {feature} match"}
    rows = [
        {
            "reason": key,
            "label": labels.get(key, key.replace("_", " ")),
            "reads": reads,
            "frac_of_reads": _div(reads, mapping["total_reads"]),
            "frac_of_unmapped": unmapped[f"{key}_frac"],
        }
        for key, reads in unmapped.items()
        if not key.endswith("_frac")
    ]
    return sorted(rows, key=lambda r: -r["reads"])


def run_summary(results: list[dict], meta: dict, cyto_outdir: str, feature: str) -> dict[str, Any]:
    """Run-level metrics every workflow reports: reads, mapping, saturation, run info."""
    recs = [r["rec"] for r in results]
    mapped = sum(r["mapped_reads"] for r in recs)
    umis = sum(r["umis"] for r in recs)
    corrected = sum(r["umi_counts"][0] for r in results)
    total_umis = sum(r["umi_counts"][1] for r in results)
    mapping = meta["mapping"]
    lib = {d["name"]: d for d in meta["library"]}
    reasons = unmapped_reasons(mapping, feature)
    return {
        "cyto_outdir": os.path.abspath(cyto_outdir),
        "total_reads": mapping["total_reads"],
        "mapped_reads": mapping["mapped_reads"],
        "mapped_reads_frac": mapping["mapped_reads_frac"],
        "top_unmapped_reason": reasons[0]["label"] if reasons else None,
        "failed_umi_qual_of_total": _div(mapping["unmapped"]["failed_umi_qual"], mapping["total_reads"]),
        "probe_barcodes_in_library": lib["probe"]["total_elem"],
        "probe_barcodes_with_reads": len(recs),
        "whitelist_size": lib["whitelist"]["total_elem"],
        "seq_saturation": 1 - umis / mapped if mapped else None,
        "umi_corrected_frac": _div(corrected, total_umis),
        "mapping_sec": meta["timings"]["Mapping"],
        "n_inputs": len(meta["run"]),
    }
