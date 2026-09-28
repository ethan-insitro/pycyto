"""Tests for ``pycyto qc`` (``pycyto.qc``).

Builds a small synthetic cyto GEX output directory at test time:

* ``A-A01``: good probe barcode, CSR filtered h5ad
* ``A-A02``: good probe barcode, CSC filtered h5ad
* ``B-H12``: no cells (no filtered h5ad; "all barcodes below rejection boundary")
* ``C-D07``: near-empty probe barcode where cyto's retainment boundary (8 UMIs) is
  below its rejection boundary (500 UMIs), i.e. an unreliable cell call
"""

import json
import os
import re

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import polars as pl
import pyarrow as pa
import pytest
import scipy.sparse as sp
from typer.testing import CliRunner

from pycyto.__main__ import app
from pycyto.qc import build_report, collect
from pycyto.qc import h5ad as qc_h5ad
from pycyto.qc.parse import parse_cyto_log, parse_filter_log, probe_sort_key

N_GENES = 40

# probe, n_cells, n_background, sparse format, near-empty
PROBES = [
    ("A-A01", 120, 1200, "csr", False),
    ("A-A02", 80, 900, "csc", False),
    ("B-H12", 0, 600, "csr", False),
    ("C-D07", 40, 500, "csr", True),
]


def _write_zst_tsv(path: str, df: pl.DataFrame) -> None:
    with pa.CompressedOutputStream(path, "zstd") as out:
        out.write(df.write_csv(separator="\t").encode())


def _write_h5ad(path: str, barcodes: np.ndarray, counts: np.ndarray, probe: str, fmt: str):
    X = sp.csr_matrix(counts) if fmt == "csr" else sp.csc_matrix(counts)
    adata = ad.AnnData(
        X=X.astype(np.float32),
        obs=pd.DataFrame(index=[f"{b}-{probe}" for b in barcodes]),
        var=pd.DataFrame(index=[f"gene{i}" for i in range(counts.shape[1])]),
    )
    adata.write_h5ad(path)


@pytest.fixture(scope="module")
def cyto_dir(tmp_path_factory):
    """Synthetic cyto output directory plus the expected per-probe metrics."""
    root = tmp_path_factory.mktemp("cyto_out")
    for sub in ("stats/reads", "stats/umi", "stats/filtering", "counts"):
        os.makedirs(root / sub, exist_ok=True)
    rng = np.random.default_rng(7)
    truth = {}
    for probe, n_cells, n_bg, fmt, near_empty in PROBES:
        n = n_cells + n_bg
        barcodes = np.array(["".join(rng.choice(list("ACGT"), 16)) for _ in range(n)])
        depth = np.r_[
            rng.lognormal(3.0 if near_empty else 8.0, 0.4, n_cells),
            rng.lognormal(2.5, 0.8, n_bg),
        ]
        counts = np.stack(
            [rng.multinomial(int(d), rng.dirichlet(np.full(N_GENES, 0.1))) for d in depth]
        )
        umis = counts.sum(1)
        reads = (umis * rng.uniform(1.1, 1.6, n)).astype(np.int64)
        _write_zst_tsv(
            str(root / "stats" / "reads" / f"{probe}.reads.tsv.zst"),
            pl.DataFrame({"barcode": barcodes, "n_umis": umis, "n_reads": reads}),
        )
        (root / "stats" / "umi" / f"{probe}.umi.json").write_text(
            json.dumps({"total": int(reads.sum()), "corrected": 5, "fraction_corrected": 5 / reads.sum()})
        )
        if n_cells:
            cells = np.argsort(-umis, kind="stable")[: n_cells - 3]
            _write_h5ad(str(root / "counts" / f"{probe}.filt.h5ad"), barcodes[cells], counts[cells], probe, fmt)
            retain = 8 if near_empty else 1500
            (root / "stats" / "filtering" / f"{probe}.log").write_text(
                f"Retainment boundary: {retain} UMIs (30 auto-accepted cells)\n"
                "Rejection boundary: 500 UMIs\n"
                "Identified 5 passing candidates and 30 retained cells.\n"
                f"Final number of filtered cells: {len(cells)}\n"
            )
            truth[probe] = {
                "cells": len(cells),
                "median_umis_per_cell": float(np.median(umis[cells])),
                "median_genes_per_cell": float(np.median((counts[cells] > 0).sum(1))),
                "total_genes_detected": int(((counts[cells] > 0).sum(0) > 0).sum()),
                "reads_in_cells": int(reads[cells].sum()),
                "mapped_reads": int(reads.sum()),
            }
        else:
            (root / "stats" / "filtering" / f"{probe}.log").write_text(
                "All barcodes have less than 500 UMIs.\nReturning empty anndata\n"
            )
            truth[probe] = {"cells": 0, "mapped_reads": int(reads.sum())}

    (root / "stats" / "mapping_map.json").write_text(
        json.dumps(
            {
                "total_reads": 1_000_000,
                "mapped_reads": 600_000,
                "unmapped_reads": 400_000,
                "mapped_reads_frac": 0.6,
                "unmapped_reads_frac": 0.4,
                "unmapped": {
                    "missing_feature": 300_000,
                    "failed_umi_qual": 120_000,
                    "missing_feature_frac": 0.75,
                    "failed_umi_qual_frac": 0.3,
                },
            }
        )
    )
    (root / "stats" / "mapping_lib.json").write_text(
        json.dumps([{"name": "probe", "total_elem": 384, "total_aggr": 384, "mate": "R1", "position": 38, "window": 5, "exact": False}])
    )
    (root / "stats" / "mapping_run.json").write_text(json.dumps([{"input_id": 0, "elapsed_sec": 12.5}]))
    (root / ".timings.tsv").write_text(
        "ibu_name\tmodule\telapsed\nAll-Barcodes\tMapping\t10.5\nA-A01\tCounting\t1.0\nA-A02\tCounting\t2.0\n"
    )
    (root / "cyto.log").write_text(
        "[2026-09-25T07:00:27.456Z INFO  cyto] Initializing...\n"
        "[2026-09-25T07:00:34.842Z INFO  cyto_workflow::gex] Running GEX Mapping Workflow\n"
        "[2026-09-25T07:00:34.842Z INFO  cyto_map::run] Using preset (GexV2) geometry: `[barcode][umi:12][:10][probe] | [gex]`\n"
        "[2026-09-25T08:11:16.440Z WARN  cyto_workflow::utils] Missing filtered h5ad file\n"
    )
    return str(root), truth


# --------------------------------------------------------------------------------------
# parsers
# --------------------------------------------------------------------------------------
class TestParsers:
    def test_filter_log_emptydrops(self, tmp_path):
        p = tmp_path / "x.log"
        p.write_text(
            "Retainment boundary: 2122 UMIs (13807 auto-accepted cells)\nRejection boundary: 500 UMIs\n"
            "Identified 778 passing candidates and 13807 retained cells.\nFinal number of filtered cells: 14585\n"
        )
        out = parse_filter_log(str(p))
        assert out["cyto_status"] == "EmptyDrops-style filtering"
        assert out["cyto_retain_umis"] == 2122
        assert out["cyto_reject_umis"] == 500
        assert out["cyto_auto_accepted"] == 13807
        assert out["cyto_passing_candidates"] == 778

    def test_filter_log_all_below(self, tmp_path):
        p = tmp_path / "x.log"
        p.write_text("All barcodes have less than 500 UMIs.\nReturning empty anndata\n")
        assert parse_filter_log(str(p))["cyto_status"] == "all barcodes below rejection boundary"

    def test_filter_log_flat_cutoff(self, tmp_path):
        p = tmp_path / "x.log"
        p.write_text(
            "Not enough barcodes to identify ambient cells. Found only 4484 barcodes.\n"
            "Returning simply filtered anndata (umis < 500)\n"
        )
        out = parse_filter_log(str(p))
        assert out["cyto_status"] == "too few barcodes; flat UMI cutoff"
        assert out["cyto_reject_umis"] == 500

    def test_filter_log_missing(self, tmp_path):
        assert parse_filter_log(str(tmp_path / "nope.log"))["cyto_status"] == "no log"

    def test_cyto_log(self, cyto_dir):
        root, _ = cyto_dir
        info = parse_cyto_log(os.path.join(root, "cyto.log"))
        assert info["preset"] == "GexV2"
        assert info["geometry"] == "[barcode][umi:12][:10][probe] | [gex]"
        assert info["workflow"] == "GEX Mapping"
        assert info["warnings"] == 1 and info["errors"] == 0
        assert info["wall_sec"] == pytest.approx(4248.984, abs=1e-3)

    def test_probe_sort_key(self):
        probes = ["B-A01", "A-B01", "A-A10", "A-A02", "BC002", "BC001"]
        assert sorted(probes, key=probe_sort_key) == ["A-A02", "A-A10", "A-B01", "B-A01", "BC001", "BC002"]


# --------------------------------------------------------------------------------------
# h5ad reader
# --------------------------------------------------------------------------------------
class TestH5ad:
    @pytest.mark.parametrize("probe", ["A-A01", "A-A02"])  # CSR and CSC
    def test_gene_stats_match_dense(self, cyto_dir, probe):
        root, _ = cyto_dir
        path = os.path.join(root, "counts", f"{probe}.filt.h5ad")
        dense = ad.read_h5ad(path).X.toarray()
        with h5py.File(path, "r") as f:
            gpc, detected = qc_h5ad.gene_stats(f)
            gpc_small, detected_small = qc_h5ad.gene_stats(f, chunk_nnz=37)  # force many chunks
            names = qc_h5ad.obs_names(f)
        np.testing.assert_array_equal(gpc, (dense > 0).sum(1))
        np.testing.assert_array_equal(detected, (dense > 0).any(0))
        np.testing.assert_array_equal(gpc, gpc_small)
        np.testing.assert_array_equal(detected, detected_small)
        assert all(n.endswith(f"-{probe}") for n in names)
        assert all(len(b) == 16 for b in qc_h5ad.strip_probe_suffix(names))


# --------------------------------------------------------------------------------------
# end to end
# --------------------------------------------------------------------------------------
class TestCollect:
    @pytest.fixture(scope="class")
    def payload(self, cyto_dir):
        root, _ = cyto_dir
        return collect(root, threads=1)

    def test_per_probe_metrics(self, payload, cyto_dir):
        _, truth = cyto_dir
        by_probe = {r["probe"]: r for r in payload["probes"]}
        assert list(by_probe) == ["A-A01", "A-A02", "B-H12", "C-D07"]
        for probe, expected in truth.items():
            rec = by_probe[probe]
            for key, value in expected.items():
                assert rec[key] == value, (probe, key)
        assert by_probe["B-H12"]["has_filtered_h5ad"] is False
        assert by_probe["B-H12"]["median_umis_per_cell"] is None

    def test_summary(self, payload, cyto_dir):
        _, truth = cyto_dir
        s = payload["summary"]
        assert s["estimated_cells"] == sum(t["cells"] for t in truth.values())
        assert s["probe_barcodes_with_cells"] == 3
        assert s["n_unreliable_probes"] == 1
        assert s["cells_in_unreliable_probes"] == truth["C-D07"]["cells"]
        assert s["mean_reads_per_cell"] == pytest.approx(1_000_000 / s["estimated_cells"])
        assert s["top_unmapped_reason"] == "no gene probe match"
        assert s["mapping_sec"] == 10.5
        assert s["total_genes_detected"] <= N_GENES

    def test_alerts(self, payload):
        titles = {a["title"] for a in payload["alerts"]}
        assert "Low fraction of reads mapped" in titles
        assert "Many reads failed UMI quality" in titles
        assert "cyto cell calling looks unreliable" in titles
        assert "Low median UMIs per cell" in titles

    def test_parallel_matches_serial(self, payload, cyto_dir):
        root, _ = cyto_dir
        parallel = collect(root, threads=2)
        assert parallel["probes"] == payload["probes"]


class TestBuildReport:
    def test_writes_html_and_csvs(self, cyto_dir, tmp_path):
        root, _ = cyto_dir
        out = str(tmp_path / "report.html")
        build_report(root, output=out, title="unit-test", threads=1)
        html = open(out).read()
        m = re.search(r'<script id="data" type="application/json">(.*?)</script>', html, re.S)
        assert m is not None
        data = json.loads(m.group(1))
        assert data["title"] == "unit-test"
        assert "__DATA__" not in html and "__TITLE__" not in html
        summary = pl.read_csv(str(tmp_path / "report_metrics_summary.csv"))
        probes = pl.read_csv(str(tmp_path / "report_probe_metrics.csv"))
        assert summary.height == 1
        assert probes.height == len(PROBES)

    def test_cli(self, cyto_dir, tmp_path):
        root, _ = cyto_dir
        out = str(tmp_path / "cli.html")
        result = CliRunner().invoke(app, ["qc", root, "--output", out, "--no-csv", "--threads", "1"])
        assert result.exit_code == 0, result.output
        assert os.path.exists(out)
        assert not os.path.exists(str(tmp_path / "cli_metrics_summary.csv"))

    def test_cli_rejects_non_cyto_dir(self, tmp_path):
        result = CliRunner().invoke(app, ["qc", str(tmp_path)])
        assert result.exit_code == 1
