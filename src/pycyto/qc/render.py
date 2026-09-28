"""Render the QC payload into the self-contained HTML report and CSV tables."""

import html
import json
import math
from importlib import resources
from typing import Any

import numpy as np
import polars as pl

TEMPLATE = "report.html"


def to_jsonable(o: Any) -> Any:
    """Make an object JSON-safe: numpy -> python, NaN/inf -> None."""
    if isinstance(o, dict):
        return {str(k): to_jsonable(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [to_jsonable(v) for v in o]
    if isinstance(o, np.ndarray):
        return to_jsonable(o.tolist())
    if isinstance(o, np.bool_):
        return bool(o)
    if isinstance(o, np.integer):
        return int(o)
    if isinstance(o, (float, np.floating)):
        f = float(o)
        return None if (math.isnan(f) or math.isinf(f)) else f
    return o


def render_html(payload: dict) -> str:
    template = resources.files("pycyto.qc").joinpath(TEMPLATE).read_text(encoding="utf-8")
    # "</" must not appear inside the <script> block that carries the data
    data = json.dumps(to_jsonable(payload), allow_nan=False).replace("</", "<\\/")
    return template.replace("__TITLE__", html.escape(payload["title"])).replace("__DATA__", data)


def write_csvs(summary: dict, probe_recs: list[dict], stem: str) -> tuple[str, str]:
    """Write ``<stem>_metrics_summary.csv`` (one row) and ``<stem>_probe_metrics.csv``."""
    summary_path = f"{stem}_metrics_summary.csv"
    probes_path = f"{stem}_probe_metrics.csv"
    pl.DataFrame([to_jsonable(summary)], infer_schema_length=None).write_csv(summary_path)
    pl.DataFrame(to_jsonable(probe_recs), infer_schema_length=None).write_csv(probes_path)
    return summary_path, probes_path
