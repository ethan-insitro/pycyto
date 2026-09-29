"""Readers for the structured outputs cyto writes next to its counts.

They assume a completed cyto run: every file is present and well-formed.
"""

import json
import os
from typing import Any

import polars as pl

from ..config import FLEX_V1_BARCODES, FLEX_V2_BARCODES

# Rank of each known barcode within its format, so `probe_sort_key` groups by prefix (V1)
# or plate position (V2). FLEX_V1_BARCODES is generated prefix-interleaved
# (BC001, CR001, AB001, BC002, …); sorting it gives the intuitive AB001..AB016,
# BC001..BC016, CR001..CR016 grouping.
_FLEX_V1_RANK = {bc: i for i, bc in enumerate(sorted(FLEX_V1_BARCODES))}
_FLEX_V2_RANK = {bc: i for i, bc in enumerate(FLEX_V2_BARCODES)}


def load_json(path: str) -> Any:
    with open(path) as fh:
        return json.load(fh)


def read_barcode_stats(path: str) -> pl.DataFrame:
    """``stats/reads/<probe>.reads.tsv.zst`` -> columns ``barcode, n_umis, n_reads``."""
    return pl.read_csv(
        path,
        separator="\t",
        schema_overrides={"barcode": pl.String, "n_umis": pl.Int64, "n_reads": pl.Int64},
    )


def load_run_metadata(cyto_outdir: str) -> dict[str, Any]:
    """Everything run-level: mapping stats and reference libraries."""
    stats = os.path.join(cyto_outdir, "stats")
    return {
        "mapping": load_json(os.path.join(stats, "mapping_map.json")),
        "library": load_json(os.path.join(stats, "mapping_lib.json")),
    }


def detect_workflow(meta: dict[str, Any]) -> str:
    """The cyto workflow (``gex``), from the feature library named in ``mapping_lib.json``."""
    libraries = {d["name"] for d in meta["library"]}
    if "gex" in libraries:
        return "gex"
    raise ValueError(
        f"Can't tell which cyto workflow produced this directory (libraries: {sorted(libraries)}); "
        "expected a `cyto workflow gex` run"
    )


def probe_sort_key(probe: str) -> tuple:
    """Flex-V2 barcodes in plate order, then Flex-V1 barcodes by prefix, then unknowns."""
    v2_rank = _FLEX_V2_RANK.get(probe.replace("_", "-"))
    if v2_rank is not None:
        return (0, v2_rank)
    v1_rank = _FLEX_V1_RANK.get(probe)
    if v1_rank is not None:
        return (1, v1_rank)
    return (2, probe)


def discover_probes(cyto_outdir: str) -> list[str]:
    reads_dir = os.path.join(cyto_outdir, "stats", "reads")
    suffix = ".reads.tsv.zst"
    names = [f.removesuffix(suffix) for f in os.listdir(reads_dir) if f.endswith(suffix)]
    return sorted(names, key=probe_sort_key)
