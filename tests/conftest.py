"""Synthetic cyto output directories for the ``pycyto qc`` tests.

The GEX run (``cyto_dir``, see ``qc_helpers.PROBES``):

* ``A-A01``: good probe barcode, CSR filtered h5ad
* ``A-A02``: good probe barcode, CSC filtered h5ad
* ``B-H12``: no cells (no filtered h5ad)
* ``C-D07``: near-empty probe barcode (cells with only ~20 UMIs)

The CRISPR run (``crispr_dir``, see ``qc_helpers.CRISPR_PROBES``) has three probe barcodes
and a guide library where some guides never appear.
"""

import json
import os

import anndata as ad
import numpy as np
import pandas as pd
import polars as pl
import pytest
import scipy.sparse as sp

from .qc_helpers import (
    CRISPR_PROBES,
    N_GENES,
    N_GUIDES,
    PROBES,
    UNUSED_GUIDES,
    write_h5ad,
    write_zst_tsv,
)


@pytest.fixture(scope="session")
def cyto_dir(tmp_path_factory):
    """Synthetic cyto output directory plus the expected per-probe metrics."""
    root = tmp_path_factory.mktemp("cyto_out")
    for sub in ("stats/reads", "stats/umi", "counts"):
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
        write_zst_tsv(
            str(root / "stats" / "reads" / f"{probe}.reads.tsv.zst"),
            pl.DataFrame({"barcode": barcodes, "n_umis": umis, "n_reads": reads}),
        )
        (root / "stats" / "umi" / f"{probe}.umi.json").write_text(
            json.dumps({"total": int(reads.sum()), "corrected": 5, "fraction_corrected": 5 / reads.sum()})
        )
        if n_cells:
            cells = np.argsort(-umis, kind="stable")[: n_cells - 3]
            write_h5ad(str(root / "counts" / f"{probe}.filt.h5ad"), barcodes[cells], counts[cells], probe, fmt)
            truth[probe] = {
                "cells": len(cells),
                "median_umis_per_cell": float(np.median(umis[cells])),
                "reads_in_cells": int(reads[cells].sum()),
                "mapped_reads": int(reads.sum()),
            }
        else:
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


@pytest.fixture(scope="session")
def crispr_dir(tmp_path_factory):
    """Synthetic ``cyto workflow crispr`` output plus the expected guide metrics."""
    root = tmp_path_factory.mktemp("crispr_out")
    for sub in ("stats/reads", "stats/umi", "counts"):
        os.makedirs(root / sub, exist_ok=True)
    rng = np.random.default_rng(11)
    guides = [f"GENE{i}_{i}" for i in range(N_GUIDES)]
    weights = np.r_[np.zeros(UNUSED_GUIDES), rng.gamma(2.0, 1.0, N_GUIDES - UNUSED_GUIDES)]
    weights[-1] *= 20  # one dominant guide
    totals = np.zeros(N_GUIDES)
    truth = {}
    for probe, n in CRISPR_PROBES:
        barcodes = np.array(["".join(rng.choice(list("ACGT"), 16)) for _ in range(n)])
        counts = np.stack([rng.multinomial(int(d), weights / weights.sum()) for d in rng.lognormal(3, 1, n)])
        umis = counts.sum(1)
        reads = (umis * rng.uniform(1.1, 1.6, n)).astype(np.int64)
        write_zst_tsv(
            str(root / "stats" / "reads" / f"{probe}.reads.tsv.zst"),
            pl.DataFrame({"barcode": barcodes, "n_umis": umis, "n_reads": reads}),
        )
        (root / "stats" / "umi" / f"{probe}.umi.json").write_text(
            json.dumps({"total": int(reads.sum()), "corrected": 3, "fraction_corrected": 3 / reads.sum()})
        )
        adata = ad.AnnData(
            X=sp.csr_matrix(counts.astype(np.float32)),
            obs=pd.DataFrame(index=[f"{b}-{probe}" for b in barcodes]),
            var=pd.DataFrame(index=guides),
        )
        adata.write_h5ad(str(root / "counts" / f"{probe}.h5ad"))
        totals += counts.sum(0)
        truth[probe] = {
            "umis": int(umis.sum()),
            "mapped_reads": int(reads.sum()),
        }
    (root / "stats" / "mapping_map.json").write_text(
        json.dumps({"total_reads": 50_000, "mapped_reads": 40_000, "mapped_reads_frac": 0.8,
                    "unmapped": {"missing_feature": 6_000, "failed_umi_qual": 4_000,
                                 "missing_feature_frac": 0.6, "failed_umi_qual_frac": 0.4}})
    )
    (root / "stats" / "mapping_lib.json").write_text(
        json.dumps([{"name": "probe", "total_elem": 384}, {"name": "whitelist", "total_elem": 737_280},
                    {"name": "crispr", "total_elem": N_GUIDES}])
    )
    (root / "stats" / "mapping_run.json").write_text(json.dumps([{"input_id": 0, "elapsed_sec": 1.5}]))
    (root / ".timings.tsv").write_text("ibu_name\tmodule\telapsed\nAll-Barcodes\tMapping\t1.2\n")
    return str(root), truth, totals, guides
