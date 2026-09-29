"""QC building blocks shared by every cyto workflow."""

import anndata as ad
import numpy as np
import polars as pl
import scipy.sparse as sp

# log10 bin edges (width 0.05) for the histograms embedded in the report
LOG_BINS = np.round(np.arange(0, 6.05, 0.05), 2)


def read_counts(path: str, chunk_rows: int = 10_000) -> tuple[pl.DataFrame, np.ndarray, list[str]]:
    """Scan a cyto count h5ad.

    Returns a ``barcode, n_features`` frame (features with nonzero counts per barcode),
    the per-feature count totals, and the feature names. ``X`` is read in row chunks
    from a backed AnnData so memory stays bounded.
    """
    adata = ad.read_h5ad(path, backed="r")
    try:
        n_features = np.zeros(adata.n_obs, dtype=np.int64)
        totals = np.zeros(adata.n_vars, dtype=np.float64)
        for start in range(0, adata.n_obs, chunk_rows):
            block = sp.csr_matrix(adata.X[start : start + chunk_rows])
            block.eliminate_zeros()
            n_features[start : start + block.shape[0]] = np.diff(block.indptr)
            totals += np.asarray(block.sum(axis=0)).ravel()
        # cyto names cells ``<barcode>-<probe>`` (e.g. ``ACGT...-A-A02``)
        barcodes = [name.split("-", 1)[0] for name in adata.obs_names]
        names = adata.var_names.tolist()
    finally:
        adata.file.close()
    return pl.DataFrame({"barcode": barcodes, "n_features": n_features}), totals, names


def rank_curve(umis_desc: np.ndarray, is_cell_desc: np.ndarray, n_points: int = 300) -> list:
    """Barcode-rank curve downsampled to ~n_points log-spaced ranks.

    Each point is ``[rank, umis, fraction of barcodes in the segment that are cells]``.
    """
    n = len(umis_desc)
    if n == 0:
        return []
    ranks = np.unique(np.clip(np.round(np.logspace(0, np.log10(n), n_points)), 1, n)).astype(int)
    prev = np.concatenate([[0], ranks[:-1]])
    csum = np.concatenate([[0], np.cumsum(is_cell_desc)])
    frac = (csum[ranks] - csum[prev]) / (ranks - prev)
    return [[int(r), int(u), round(float(f), 3)] for r, u, f in zip(ranks, umis_desc[ranks - 1], frac)]


def log_hist(values: np.ndarray) -> list[int] | None:
    if len(values) == 0:
        return None
    counts, _ = np.histogram(np.log10(np.maximum(values, 1)), bins=np.append(LOG_BINS, 6.05))
    return counts.tolist()
