"""``pycyto qc`` on a CRISPR run."""

import numpy as np
import pytest

from pycyto.qc import build_report, collect

from .qc_helpers import N_GUIDES, UNUSED_GUIDES


@pytest.fixture(scope="module")
def payload(crispr_dir):
    return collect(crispr_dir[0])


class TestCrispr:
    def test_per_probe_metrics(self, payload, crispr_dir):
        _, truth, _, _ = crispr_dir
        by_probe = {r["probe"]: r for r in payload["probes"]}
        assert payload["workflow"] == "crispr"
        for probe, expected in truth.items():
            for key, value in expected.items():
                assert by_probe[probe][key] == value, (probe, key)
        assert all(r["flag"] is None for r in payload["probes"])

    def test_summary(self, payload, crispr_dir):
        _, _, totals, _ = crispr_dir
        s = payload["summary"]
        assert s["guides_in_library"] == N_GUIDES
        assert s["guides_detected"] == int((totals > 0).sum()) == N_GUIDES - UNUSED_GUIDES
        assert s["frac_guides_detected"] == pytest.approx((N_GUIDES - UNUSED_GUIDES) / N_GUIDES)
        assert s["guide_umis"] == int(totals.sum())
        assert s["median_umis_per_guide"] == float(np.median(totals))
        p10, p90 = np.percentile(totals, [10, 90])
        assert s["guide_skew_ratio"] == (p90 / p10 if p10 else None)
        assert s["top_unmapped_reason"] == "No guide match"

    def test_top_guides(self, payload, crispr_dir):
        _, _, totals, guides = crispr_dir
        top = payload["pooled"]["top_guides"]
        assert top[0]["guide"] == guides[int(np.argmax(totals))]
        assert top[0]["umis"] == int(totals.max())
        assert [g["umis"] for g in top] == sorted((g["umis"] for g in top), reverse=True)

    def test_alerts(self, payload):
        titles = {a["title"] for a in payload["alerts"]}
        assert "Guides missing from the library" in titles  # 80% < 90% expected
        assert "Low fraction of reads in cells" not in titles  # GEX-only alert

    def test_report(self, crispr_dir, tmp_path):
        out = str(tmp_path / "crispr.html")
        build_report(crispr_dir[0], output=out)
        html = open(out).read()
        assert "subtitle: \"cyto workflow crispr" in html and "Estimated number of cells" not in html
        assert html.rstrip().endswith("</html>") and "<script>main();</script>" in html
