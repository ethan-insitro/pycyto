"""Parsers for the metadata files cyto writes next to its counts.

All parsers are tolerant: a missing or malformed file yields an empty/default value
so that a partial cyto output directory still produces a report.
"""

import datetime as dt
import json
import logging
import os
import re
from typing import Any

import numpy as np
import polars as pl

from ..config import _is_flex_v2_barcode

logger = logging.getLogger("pycyto.qc")

# unmapped-read categories from ``stats/mapping_map.json`` -> human-readable labels
UNMAPPED_LABELS = {
    "missing_feature": "no gene probe match",
    "missing_probe": "no probe barcode match",
    "failed_umi_qual": "UMI failed quality",
    "missing_whitelist": "cell barcode not in whitelist",
    "umi_truncated": "UMI truncated",
}


def load_json(path: str, default: Any = None) -> Any:
    try:
        with open(path) as fh:
            return json.load(fh)
    except (OSError, json.JSONDecodeError):
        logger.debug(f"Could not read JSON: {path}")
        return default


def read_barcode_stats(path: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Read ``stats/reads/<probe>.reads.tsv.zst`` -> (barcodes, n_umis, n_reads)."""
    df = pl.read_csv(
        path,
        separator="\t",
        has_header=True,
        schema_overrides={"barcode": pl.String, "n_umis": pl.Int64, "n_reads": pl.Int64},
    )
    return (
        df["barcode"].to_numpy().astype(str),
        df["n_umis"].to_numpy().astype(np.int64),
        df["n_reads"].to_numpy().astype(np.int64),
    )


def parse_filter_log(path: str) -> dict:
    """Parse cyto's droplet-filtering log (``stats/filtering/<probe>.log``)."""
    out: dict[str, Any] = {
        "cyto_status": "no log",
        "cyto_retain_umis": None,
        "cyto_reject_umis": None,
        "cyto_auto_accepted": None,
        "cyto_passing_candidates": None,
    }
    try:
        with open(path) as fh:
            txt = fh.read()
    except OSError:
        return out
    if "Final number of filtered cells" in txt:
        out["cyto_status"] = "EmptyDrops-style filtering"
    if "All barcodes have less than" in txt:
        out["cyto_status"] = "all barcodes below rejection boundary"
    elif "Not enough barcodes to identify ambient" in txt:
        out["cyto_status"] = "too few barcodes; flat UMI cutoff"
        if m := re.search(r"Returning simply filtered anndata \(umis < (\d+)\)", txt):
            out["cyto_reject_umis"] = int(m.group(1))
    if m := re.search(r"Retainment boundary:\s*(\d+) UMIs \((\d+) auto-accepted", txt):
        out["cyto_retain_umis"] = int(m.group(1))
        out["cyto_auto_accepted"] = int(m.group(2))
    if m := re.search(r"Rejection boundary:\s*(\d+) UMIs", txt):
        out["cyto_reject_umis"] = int(m.group(1))
    if m := re.search(r"Identified (\d+) passing candidates", txt):
        out["cyto_passing_candidates"] = int(m.group(1))
    return out


def parse_cyto_log(path: str) -> dict:
    """Pull geometry, workflow, timestamps and warning/error counts from ``cyto.log``."""
    info: dict[str, Any] = {
        "geometry": None,
        "preset": None,
        "workflow": None,
        "start": None,
        "end": None,
        "wall_sec": None,
        "warnings": 0,
        "errors": 0,
    }
    try:
        with open(path) as fh:
            lines = fh.read().splitlines()
    except OSError:
        return info
    ts_re = re.compile(r"^\[(\d{4}-\d\d-\d\dT[\d:.]+)Z\s+(\w+)")
    stamps: list[dt.datetime] = []
    for ln in lines:
        if m := ts_re.match(ln):
            try:
                stamps.append(dt.datetime.fromisoformat(m.group(1)[:26]))
            except ValueError:
                pass
            if m.group(2) == "WARN":
                info["warnings"] += 1
            elif m.group(2) == "ERROR":
                info["errors"] += 1
        if g := re.search(r"Using (?:preset \((\w+)\) )?geometry: `([^`]+)`", ln):
            info["preset"], info["geometry"] = g.group(1), g.group(2)
        if w := re.search(r"Running (.+?) Workflow", ln):
            info["workflow"] = w.group(1)
    if stamps:
        first, last = min(stamps), max(stamps)
        info["start"] = first.isoformat(sep=" ", timespec="seconds")
        info["end"] = last.isoformat(sep=" ", timespec="seconds")
        info["wall_sec"] = (last - first).total_seconds()
    return info


def parse_timings(path: str) -> dict[str, float]:
    """Sum ``.timings.tsv`` elapsed seconds by pipeline module."""
    if not os.path.exists(path):
        return {}
    try:
        df = pl.read_csv(path, separator="\t")
    except Exception:  # malformed file shouldn't sink the report
        logger.debug(f"Could not parse timings: {path}")
        return {}
    totals: dict[str, float] = {}
    for module, elapsed in df.select(["module", "elapsed"]).iter_rows():
        totals[str(module)] = totals.get(str(module), 0.0) + float(elapsed)
    return totals


def probe_sort_key(probe: str) -> tuple:
    """Sort Flex-V2 barcodes by set/row/column and everything else lexically."""
    if _is_flex_v2_barcode(probe):
        return (0, probe[0], probe[2], int(probe[3:]))
    return (1, probe, "", 0)


def discover_probes(cyto_outdir: str) -> list[str]:
    reads_dir = os.path.join(cyto_outdir, "stats", "reads")
    suffix = ".reads.tsv.zst"
    try:
        names = [f[: -len(suffix)] for f in os.listdir(reads_dir) if f.endswith(suffix)]
    except OSError:
        return []
    return sorted(names, key=probe_sort_key)
