"""Write the QC metrics to CSV."""

import polars as pl


def write_csvs(payload: dict, stem: str) -> list[str]:
    """``<stem>_metrics_summary.csv`` (one row of run-level metrics)."""
    path = f"{stem}_metrics_summary.csv"
    pl.DataFrame([payload["summary"]], infer_schema_length=None).write_csv(path)
    return [path]
