"""Synthetic cyto output directories for the ``pycyto qc`` tests.

The GEX run (``cyto_dir``, see ``qc_helpers.PROBES``):

* ``A-A01``, ``A-A02``: probe barcodes with many high-UMI barcodes
* ``B-H12``: background barcodes only
* ``C-D07``: near-empty probe barcode (high-UMI barcodes with only ~20 UMIs)
"""

import json
import os

import numpy as np
import polars as pl
import pytest

from .qc_helpers import N_GENES, PROBES, write_zst_tsv


@pytest.fixture(scope="session")
def cyto_dir(tmp_path_factory):
    """Synthetic cyto output directory plus the expected per-probe metrics."""
    root = tmp_path_factory.mktemp("cyto_out")
    for sub in ("stats/reads", "stats/umi"):
        os.makedirs(root / sub, exist_ok=True)
    rng = np.random.default_rng(7)
    truth = {}
    for probe, n_cells, n_bg, _, near_empty in PROBES:
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
        write_zst_tsv(
            str(root / "stats" / "reads" / f"{probe}.reads.tsv.zst"),
            pl.DataFrame({"barcode": barcodes, "n_umis": umis, "n_reads": reads}),
        )
        (root / "stats" / "umi" / f"{probe}.umi.json").write_text(
            json.dumps({"total": int(reads.sum()), "corrected": 5, "fraction_corrected": 5 / reads.sum()})
        )
        truth[probe] = {"mapped_reads": int(reads.sum()), "umis": int(umis.sum())}

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
        json.dumps([
            {"name": "probe", "total_elem": 384, "total_aggr": 384, "mate": "R1", "position": 38, "window": 5, "exact": False},
            {"name": "whitelist", "total_elem": 737_280, "total_aggr": 737_280, "mate": "R1", "position": 0, "window": 5, "exact": False},
            {"name": "gex", "total_elem": 3 * N_GENES, "total_aggr": N_GENES, "mate": "R2", "position": 0, "window": 5, "exact": False},
        ])
    )
    (root / "stats" / "mapping_run.json").write_text(json.dumps([{"input_id": 0, "elapsed_sec": 12.5}]))
    (root / ".timings.tsv").write_text(
        "ibu_name\tmodule\telapsed\nAll-Barcodes\tMapping\t10.5\nA-A01\tCounting\t1.0\nA-A02\tCounting\t2.0\n"
    )
    return str(root), truth
