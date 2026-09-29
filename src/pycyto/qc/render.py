"""Render the QC payload into the self-contained HTML report and CSV tables."""

import html
import json
from importlib import resources

import polars as pl


def _asset(name: str) -> str:
    return resources.files("pycyto.qc").joinpath(name).read_text(encoding="utf-8")


def render_html(payload: dict) -> str:
    """Inline the CSS, the shared ``report.js``, the workflow's ``report_<workflow>.js``
    and the payload into ``report.html``."""
    # the payload lives inside a <script> block, so "</" must be escaped
    data = json.dumps(payload, allow_nan=False).replace("</", "<\\/")
    # data goes last so nothing inside the payload is mistaken for a placeholder
    return (
        _asset("report.html")
        .replace("/*__CSS__*/", _asset("report.css"))
        .replace("/*__JS__*/", _asset("report.js"))
        .replace("/*__WORKFLOW_JS__*/", _asset(f"report_{payload['workflow']}.js"))
        .replace("__TITLE__", html.escape(payload["title"]))
        .replace("__DATA__", data)
    )


def write_csvs(payload: dict, stem: str) -> list[str]:
    """``<stem>_metrics_summary.csv`` (one row of run-level metrics)."""
    path = f"{stem}_metrics_summary.csv"
    pl.DataFrame([payload["summary"]], infer_schema_length=None).write_csv(path)
    return [path]
