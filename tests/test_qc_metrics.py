"""``pycyto.qc.metrics``: shared building blocks."""

import os

import anndata as ad
import numpy as np
import pytest

from pycyto.qc.metrics import rank_curve, read_counts


class TestReadCounts:
    @pytest.mark.parametrize("probe", ["A-A01", "A-A02"])  # CSR and CSC
    def test_matches_dense(self, cyto_dir, probe):
        root, _ = cyto_dir
        path = os.path.join(root, "counts", f"{probe}.filt.h5ad")
        adata = ad.read_h5ad(path)
        dense = adata.X.toarray()
        cells, totals, names = read_counts(path)
        cells_small, totals_small, _ = read_counts(path, chunk_rows=7)  # force many chunks
        np.testing.assert_array_equal(cells["n_features"].to_numpy(), (dense > 0).sum(1))
        np.testing.assert_allclose(totals, dense.sum(0))
        assert cells.equals(cells_small)
        np.testing.assert_allclose(totals, totals_small)
        assert names == adata.var_names.tolist()
        # probe suffix stripped from obs names
        assert cells["barcode"].to_list() == [n.split("-", 1)[0] for n in adata.obs_names]
        assert cells["barcode"].str.len_chars().eq(16).all()


def test_rank_curve():
    rng = np.random.default_rng(0)
    umis = np.sort(rng.integers(1, 10_000, 5_000))[::-1]
    is_cell = umis > 2_000
    curve = rank_curve(umis, is_cell, n_points=50)
    ranks = [r for r, _, _ in curve]
    assert ranks[0] == 1 and ranks[-1] == len(umis) and ranks == sorted(set(ranks))
    prev = 0
    for rank, u, frac in curve:  # each point summarizes barcodes (prev, rank]
        assert u == umis[rank - 1]
        assert frac == round(float(is_cell[prev:rank].mean()), 3)
        prev = rank
    assert rank_curve(np.array([], dtype=int), np.array([], dtype=bool)) == []
