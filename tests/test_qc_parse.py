"""``pycyto.qc.parse``: cyto metadata readers and workflow detection."""

import pytest

from pycyto.qc.parse import (
    detect_workflow,
    load_run_metadata,
    probe_sort_key,
)


class TestParsers:
    def test_probe_sort_key(self):
        probes = ["B-A01", "A-B01", "A-A10", "A-A02", "BC002", "BC001"]
        assert sorted(probes, key=probe_sort_key) == ["A-A02", "A-A10", "A-B01", "B-A01", "BC001", "BC002"]


def test_detect_workflow(cyto_dir):
    assert detect_workflow(load_run_metadata(cyto_dir[0])) == "gex"
    with pytest.raises(ValueError):
        detect_workflow({"library": [{"name": "probe"}]})
