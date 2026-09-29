"""``pycyto.qc.collect`` on a GEX run."""

import pytest

from pycyto.qc import collect

from .qc_helpers import N_GENES


@pytest.fixture(scope="module")
def payload(cyto_dir):
    root, _ = cyto_dir
    return collect(root)


class TestCollect:
    def test_per_probe_metrics(self, payload, cyto_dir):
        _, truth = cyto_dir
        by_probe = {r["probe"]: r for r in payload["probes"]}
        assert list(by_probe) == ["A-A01", "A-A02", "B-H12", "C-D07"]
        for probe, expected in truth.items():
            rec = by_probe[probe]
            for key, value in expected.items():
                assert rec[key] == value, (probe, key)
        assert by_probe["B-H12"]["median_umis_per_cell"] is None

    def test_summary(self, payload, cyto_dir):
        _, truth = cyto_dir
        s = payload["summary"]
        assert s["estimated_cells"] == sum(t["cells"] for t in truth.values())
        assert s["probe_barcodes_with_cells"] == 3
        assert s["mean_reads_per_cell"] == pytest.approx(1_000_000 / s["estimated_cells"])
        assert s["top_unmapped_reason"] == "No gene probe match"
        assert s["total_genes_detected"] <= N_GENES

    def test_alerts(self, payload):
        titles = {a["title"] for a in payload["alerts"]}
        assert "Low fraction of reads mapped" in titles
        assert "Many reads failed UMI quality" in titles
        assert "Low median UMIs per cell" in titles
