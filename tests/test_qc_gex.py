"""``pycyto.qc.collect`` on a GEX run."""

import pytest

from pycyto.qc import collect


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

    def test_summary(self, payload, cyto_dir):
        _, truth = cyto_dir
        s = payload["summary"]
        assert s["total_reads"] == 1_000_000
        assert s["mapped_reads_frac"] == 0.6
        assert s["probe_barcodes_with_reads"] == 4
        assert s["probe_barcodes_in_library"] == 384
        umis = sum(t["umis"] for t in truth.values())
        mapped = sum(t["mapped_reads"] for t in truth.values())
        assert s["seq_saturation"] == pytest.approx(1 - umis / mapped)
