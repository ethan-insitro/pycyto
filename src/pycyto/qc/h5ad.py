"""Minimal, memory-bounded readers for cyto's filtered h5ad files.

The QC only needs cell names and per-cell / per-gene detection counts, so we read
the on-disk anndata layout with h5py directly instead of materializing an AnnData
(filtered files can be several hundred MB each).
"""

import h5py
import numpy as np


def _read_str_array(node) -> np.ndarray:
    """Read an anndata string array (plain dataset or nullable-string-array group)."""
    if isinstance(node, h5py.Group):
        node = node["values"]
    arr = node[()]
    if arr.dtype.kind in ("S", "O"):
        arr = np.array([x.decode() if isinstance(x, bytes) else str(x) for x in arr])
    return arr.astype(str)


def obs_names(f) -> np.ndarray:
    obs = f["obs"]
    idx = obs.attrs.get("_index", "_index")
    if isinstance(idx, bytes):
        idx = idx.decode()
    return _read_str_array(obs[idx])


def n_vars(f) -> int:
    X = f["X"]
    if isinstance(X, h5py.Group):
        return int(X.attrs["shape"][1])
    return int(X.shape[1])


def strip_probe_suffix(names: np.ndarray) -> np.ndarray:
    """cyto names cells ``<barcode>-<probe>`` (e.g. ``ACGT...-A-A02``); keep ``<barcode>``."""
    return np.array([n.split("-", 1)[0] for n in names])


def gene_stats(f, chunk_nnz: int = 50_000_000) -> tuple[np.ndarray, np.ndarray]:
    """Genes detected per cell, and whether each gene is detected in any cell.

    Handles CSR, CSC and dense ``X``; sparse matrices are streamed in chunks of
    ``chunk_nnz`` stored values so memory stays bounded.
    """
    X = f["X"]
    nv = n_vars(f)
    if isinstance(X, h5py.Group):
        enc = X.attrs.get("encoding-type", b"csr_matrix")
        enc = enc.decode() if isinstance(enc, bytes) else str(enc)
        indptr = X["indptr"][()].astype(np.int64)
        indices, data = X["indices"], X["data"]
        n_obs = int(X.attrs["shape"][0])
        csr = enc.startswith("csr")
        n_major = n_obs if csr else nv
        per_obs = np.zeros(n_obs, dtype=np.int64)
        per_var = np.zeros(nv, dtype=np.int64)
        nnz_major = np.diff(indptr)
        i0 = 0
        while i0 < n_major:
            i1 = int(np.searchsorted(indptr, indptr[i0] + chunk_nnz, side="right")) - 1
            i1 = max(min(i1, n_major), i0 + 1)
            a, b = indptr[i0], indptr[i1]
            if b > a:
                minor = np.asarray(indices[a:b])
                major = np.repeat(np.arange(i0, i1), nnz_major[i0:i1])
                keep = np.asarray(data[a:b]) != 0
                rows, cols = (major, minor) if csr else (minor, major)
                per_obs += np.bincount(rows[keep], minlength=n_obs)
                per_var += np.bincount(cols[keep], minlength=nv)
            i0 = i1
        return per_obs, per_var > 0

    # dense
    n_obs = X.shape[0]
    gpc: list[np.ndarray] = []
    detected = np.zeros(nv, dtype=bool)
    step = max(1, chunk_nnz // max(nv, 1))
    for r0 in range(0, n_obs, step):
        blk = np.asarray(X[r0 : r0 + step]) != 0
        gpc.append(blk.sum(1))
        detected |= blk.any(0)
    return (np.concatenate(gpc) if gpc else np.zeros(0, dtype=np.int64)), detected
