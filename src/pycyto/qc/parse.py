"""Readers for the structured outputs cyto writes next to its counts.

They assume a completed cyto run: every file is present and well-formed.
"""

import json
import os
from typing import Any

import polars as pl

from ..config import _is_flex_v2_barcode


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
    for workflow in ("gex",):
        if workflow in libraries:
            return workflow
    raise ValueError(
        f"Can't tell which cyto workflow produced this directory (libraries: {sorted(libraries)}); "
        "expected a `cyto workflow gex` run"
    )


def probe_sort_key(probe: str) -> tuple:
    """Flex-V2 barcodes by set/row/column, then everything else lexically."""
    if _is_flex_v2_barcode(probe):
        return (0, probe[0], probe[2], int(probe[3:]))
    return (1, probe, "", 0)


def discover_probes(cyto_outdir: str) -> list[str]:
    reads_dir = os.path.join(cyto_outdir, "stats", "reads")
    suffix = ".reads.tsv.zst"
    names = [f.removesuffix(suffix) for f in os.listdir(reads_dir) if f.endswith(suffix)]
    return sorted(names, key=probe_sort_key)
