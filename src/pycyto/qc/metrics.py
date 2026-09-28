"""Per-probe-barcode and run-level QC metrics.

Cells are exactly the barcodes in cyto's ``counts/<probe>.filt.h5ad``; probe barcodes
without that file have no cells.
"""

import logging
import os
from typing import Any

import h5py
import numpy as np

from . import h5ad
from .parse import UNMAPPED_LABELS, load_json, parse_filter_log, read_barcode_stats

logger = logging.getLogger("pycyto.qc")

# log10 bins (width 0.05) for the UMI / gene histograms embedded in the report
LOG_BINS = np.round(np.arange(0, 6.05, 0.05), 2)


def rank_curve(
    umis_desc: np.ndarray, is_cell_desc: np.ndarray, n_points: int = 300
) -> list[list[float]]:
    """Downsample a barcode-rank curve to ~n_points log-spaced ranks.

    Each point is ``[rank, umis, fraction_of_barcodes_in_segment_that_are_cells]``.
    """
    n = len(umis_desc)
    if n == 0:
        return []
    ranks = np.unique(np.round(np.logspace(0, np.log10(n), n_points)).astype(int))
    ranks = np.clip(ranks, 1, n)
    csum = np.concatenate([[0], np.cumsum(is_cell_desc)])
    out, prev = [], 0
    for r in ranks:
        seg = r - prev
        frac = (csum[r] - csum[prev]) / seg if seg > 0 else float(is_cell_desc[r - 1])
        out.append([int(r), int(umis_desc[r - 1]), round(float(frac), 3)])
        prev = r
    return out


def log_hist(values: np.ndarray) -> list[int]:
    v = np.log10(np.maximum(values, 1))
    h, _ = np.histogram(v, bins=np.append(LOG_BINS, LOG_BINS[-1] + 0.05))
    return h.astype(int).tolist()


def process_probe(cyto_outdir: str, probe: str) -> dict[str, Any]:
    """Compute metrics for one probe barcode.

    Returns ``{"rec": <flat metrics dict>, "extras": <plot data + per-cell arrays>}``.
    """
    rec: dict[str, Any] = {"probe": probe}
    bcs, umis, reads = read_barcode_stats(
        os.path.join(cyto_outdir, "stats", "reads", f"{probe}.reads.tsv.zst")
    )
    rec["n_barcodes"] = int(len(bcs))
    rec["mapped_reads"] = int(reads.sum())
    rec["umis"] = int(umis.sum())
    rec["seq_saturation"] = float(1 - umis.sum() / reads.sum()) if reads.sum() else None

    umi_json = load_json(os.path.join(cyto_outdir, "stats", "umi", f"{probe}.umi.json"), {}) or {}
    rec["umi_corrected"] = umi_json.get("corrected")
    rec["umi_total"] = umi_json.get("total")
    rec["umi_corrected_frac"] = umi_json.get("fraction_corrected")

    rec.update(parse_filter_log(os.path.join(cyto_outdir, "stats", "filtering", f"{probe}.log")))

    filt = os.path.join(cyto_outdir, "counts", f"{probe}.filt.h5ad")
    rec["has_filtered_h5ad"] = os.path.exists(filt)
    cell_mask = np.zeros(len(bcs), dtype=bool)
    genes_per_cell = np.zeros(0, dtype=np.int64)
    detected: np.ndarray | None = None
    n_vars: int | None = None
    if rec["has_filtered_h5ad"]:
        with h5py.File(filt, "r") as f:
            cell_bcs = h5ad.strip_probe_suffix(h5ad.obs_names(f))
            gpc_h5, detected = h5ad.gene_stats(f)
            n_vars = h5ad.n_vars(f)
        order = {b: i for i, b in enumerate(cell_bcs)}
        cell_mask = np.fromiter((b in order for b in bcs), dtype=bool, count=len(bcs))
        missing = len(cell_bcs) - int(cell_mask.sum())
        if missing:
            logger.warning(f"[{probe}] - {missing} filtered barcodes missing from reads stats")
        # align genes-per-cell to the reads-stats barcode order
        idx = np.array([order[b] for b in bcs[cell_mask]], dtype=np.int64)
        genes_per_cell = gpc_h5[idx]

    n_cells = int(cell_mask.sum())
    cell_umis = umis[cell_mask]
    rec["cells"] = n_cells
    rec["reads_in_cells"] = int(reads[cell_mask].sum())
    rec["frac_reads_in_cells"] = (
        rec["reads_in_cells"] / rec["mapped_reads"] if rec["mapped_reads"] else None
    )
    rec["frac_umis_in_cells"] = float(cell_umis.sum() / umis.sum()) if umis.sum() else None
    rec["mean_reads_per_cell"] = rec["mapped_reads"] / n_cells if n_cells else None
    rec["median_umis_per_cell"] = float(np.median(cell_umis)) if n_cells else None
    rec["median_genes_per_cell"] = float(np.median(genes_per_cell)) if n_cells else None
    rec["total_genes_detected"] = int(detected.sum()) if detected is not None else None
    rec["n_vars"] = n_vars

    order_desc = np.argsort(-umis, kind="stable")
    extras = {
        "curve": rank_curve(umis[order_desc], cell_mask[order_desc]),
        "umi_hist": log_hist(cell_umis) if n_cells else None,
        "gene_hist": log_hist(genes_per_cell) if n_cells else None,
        "cell_umis": cell_umis.astype(np.int64),
        "cell_genes": genes_per_cell.astype(np.int64),
        "detected": detected,
    }
    return {"rec": rec, "extras": extras}


def _unreliable(rec: dict) -> bool:
    """cyto auto-accepted barcodes below its own rejection boundary."""
    retain, reject = rec.get("cyto_retain_umis"), rec.get("cyto_reject_umis")
    return rec["cells"] > 0 and retain is not None and reject is not None and retain < reject


def summarize(
    recs: list[dict],
    cell_umis: np.ndarray,
    cell_genes: np.ndarray,
    detected: list[np.ndarray],
    mapping: dict,
    lib: list[dict],
    loginfo: dict,
    timings: dict[str, float],
    n_inputs: int,
    cyto_outdir: str,
) -> dict[str, Any]:
    """Run-level metrics (the numbers behind the report's Summary tab)."""
    total_reads = mapping.get("total_reads")
    mapped_reads = mapping.get("mapped_reads")
    unmapped = mapping.get("unmapped", {}) or {}
    probe_mapped = sum(r["mapped_reads"] for r in recs)
    n_cells = sum(r["cells"] for r in recs)
    with_cells = [r for r in recs if r["cells"] > 0]
    unreliable = [r for r in recs if _unreliable(r)]
    umi_total = sum(r["umi_total"] or 0 for r in recs)
    umi_corrected = sum(r["umi_corrected"] or 0 for r in recs)

    reasons = {k: v for k, v in unmapped.items() if not k.endswith("_frac")}
    top_reason = max(reasons, key=lambda k: reasons[k]) if reasons else None

    lib_by = {d.get("name"): d for d in lib}
    n_vars = next((r["n_vars"] for r in recs if r.get("n_vars")), None)
    total_genes = None
    if detected and len({len(d) for d in detected}) == 1:
        total_genes = int(np.logical_or.reduce(detected).sum())

    failed_umi = unmapped.get("failed_umi_qual")
    return {
        "cyto_outdir": os.path.abspath(cyto_outdir),
        "total_reads": total_reads,
        "mapped_reads": mapped_reads,
        "mapped_reads_frac": mapping.get("mapped_reads_frac"),
        "top_unmapped_reason": UNMAPPED_LABELS.get(top_reason, str(top_reason).replace("_", " "))
        if top_reason
        else None,
        "failed_umi_qual_of_total": failed_umi / total_reads
        if (total_reads and failed_umi is not None)
        else None,
        "estimated_cells": n_cells,
        "cells_in_unreliable_probes": sum(r["cells"] for r in unreliable),
        "n_unreliable_probes": len(unreliable),
        "probe_barcodes_in_library": lib_by.get("probe", {}).get("total_elem"),
        "probe_barcodes_with_reads": len(recs),
        "probe_barcodes_with_cells": len(with_cells),
        "n_probes_without_cells": len(recs) - len(with_cells),
        "mean_reads_per_cell": total_reads / n_cells if (total_reads and n_cells) else None,
        "mean_mapped_reads_per_cell": mapped_reads / n_cells if (mapped_reads and n_cells) else None,
        "median_umis_per_cell": float(np.median(cell_umis)) if len(cell_umis) else None,
        "median_genes_per_cell": float(np.median(cell_genes)) if len(cell_genes) else None,
        "total_genes_detected": total_genes,
        "genes_in_reference": n_vars or lib_by.get("gex", {}).get("total_aggr"),
        "gex_probes_in_reference": lib_by.get("gex", {}).get("total_elem"),
        "whitelist_size": lib_by.get("whitelist", {}).get("total_elem"),
        "frac_reads_in_cells": sum(r["reads_in_cells"] for r in recs) / probe_mapped
        if probe_mapped
        else None,
        "background_probe_read_frac": sum(r["mapped_reads"] for r in recs if r["cells"] == 0)
        / probe_mapped
        if probe_mapped
        else None,
        "seq_saturation": 1 - sum(r["umis"] for r in recs) / probe_mapped if probe_mapped else None,
        "umi_corrected_frac": umi_corrected / umi_total if umi_total else None,
        "cells_median_per_probe": float(np.median([r["cells"] for r in with_cells]))
        if with_cells
        else None,
        "geometry": loginfo.get("geometry"),
        "preset": loginfo.get("preset"),
        "workflow": loginfo.get("workflow"),
        "run_start": loginfo.get("start"),
        "run_end": loginfo.get("end"),
        "run_wall_sec": loginfo.get("wall_sec"),
        "log_warnings": loginfo.get("warnings"),
        "log_errors": loginfo.get("errors"),
        "mapping_sec": timings.get("Mapping"),
        "n_inputs": n_inputs,
    }
