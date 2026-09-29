"""Per-probe and run-level QC metrics for a ``cyto workflow gex`` run."""

import os
from typing import Any

from .parse import load_json, read_barcode_stats


def _div(a, b) -> float | None:
    return a / b if a is not None and b else None


def process_probe(cyto_outdir: str, probe: str) -> dict[str, Any]:
    """Metrics for one probe barcode.

    ``rec`` is the row shown in the report's probe table; ``umi_counts`` is the
    ``(corrected, total)`` pair from ``stats/umi/<probe>.umi.json``, only used by :func:`summarize`.
    """
    stats = os.path.join(cyto_outdir, "stats")
    df = read_barcode_stats(os.path.join(stats, "reads", f"{probe}.reads.tsv.zst"))
    umi = load_json(os.path.join(stats, "umi", f"{probe}.umi.json"))
    return {
        "rec": {"probe": probe, "mapped_reads": df["n_reads"].sum(), "umis": df["n_umis"].sum()},
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
    return {
        "cyto_outdir": os.path.abspath(cyto_outdir),
        "total_reads": mapping["total_reads"],
        "mapped_reads": mapping["mapped_reads"],
        "mapped_reads_frac": mapping["mapped_reads_frac"],
        "probe_barcodes_in_library": lib["probe"]["total_elem"],
        "probe_barcodes_with_reads": len(recs),
        "seq_saturation": 1 - umis / mapped if mapped else None,
        "umi_corrected_frac": _div(corrected, total_umis),
    }
