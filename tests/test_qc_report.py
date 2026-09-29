"""``pycyto qc``: the CLI and metric CSVs."""

import os

import polars as pl
from typer.testing import CliRunner

from pycyto.__main__ import app
from pycyto.qc import build_report


class TestBuildReport:
    def test_writes_csvs(self, cyto_dir, tmp_path):
        root, _ = cyto_dir
        build_report(root, output=str(tmp_path / "report"), threads=1)
        summary = pl.read_csv(str(tmp_path / "report_metrics_summary.csv"))
        assert summary.height == 1

    def test_cli(self, cyto_dir, tmp_path):
        root, _ = cyto_dir
        out = str(tmp_path / "cli")
        result = CliRunner().invoke(app, ["qc", root, "--output", out, "--threads", "1"])
        assert result.exit_code == 0, result.output
        assert os.path.exists(f"{out}_metrics_summary.csv")

    def test_cli_rejects_non_cyto_dir(self, tmp_path):
        result = CliRunner().invoke(app, ["qc", str(tmp_path)])
        assert result.exit_code == 1
