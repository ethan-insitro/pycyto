"""``pycyto qc``: the CLI, metric CSVs and HTML report."""

import json
import os
import re

import polars as pl
from typer.testing import CliRunner

from pycyto.__main__ import app
from pycyto.qc import build_report


class TestBuildReport:
    def test_writes_html(self, cyto_dir, tmp_path):
        root, _ = cyto_dir
        out = str(tmp_path / "report.html")
        build_report(root, output=out, title="unit-test", write_csv=False)
        html = open(out).read()
        m = re.search(r'<script id="data" type="application/json">(.*?)</script>', html, re.S)
        assert m is not None
        data = json.loads(m.group(1))
        assert data["title"] == "unit-test"
        for placeholder in ("__DATA__", "__TITLE__", "/*__CSS__*/", "/*__JS__*/"):
            assert placeholder not in html
        assert "--cell:" in html and "function rankCurve" in html  # CSS and JS inlined
        assert "Estimated number of cells" in html  # the GEX workflow config, not CRISPR's
        assert data["workflow"] == "gex"
        # chart libraries come from a pinned, integrity-checked CDN URL
        for lib in ("d3@7.9.0/dist/d3.min.js", "plot@0.6.17/dist/plot.umd.min.js"):
            assert re.search(rf'<script src="[^"]*{re.escape(lib)}"\s+integrity="sha256-', html)
        assert not os.path.exists(str(tmp_path / "report_metrics_summary.csv"))

    def test_writes_csvs(self, cyto_dir, tmp_path):
        root, _ = cyto_dir
        build_report(root, output=str(tmp_path / "report.html"))
        summary = pl.read_csv(str(tmp_path / "report_metrics_summary.csv"))
        assert summary.height == 1

    def test_cli(self, cyto_dir, tmp_path):
        root, _ = cyto_dir
        out = str(tmp_path / "cli.html")
        result = CliRunner().invoke(app, ["qc", root, "--output", out, "--no-csv"])
        assert result.exit_code == 0, result.output
        assert os.path.exists(out)
        assert not os.path.exists(str(tmp_path / "cli_metrics_summary.csv"))

    def test_cli_rejects_non_cyto_dir(self, tmp_path):
        result = CliRunner().invoke(app, ["qc", str(tmp_path)])
        assert result.exit_code == 1
