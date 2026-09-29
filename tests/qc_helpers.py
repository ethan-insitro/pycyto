"""Constants and writers for the synthetic cyto runs in ``conftest.py``."""

import anndata as ad
import numpy as np
import pandas as pd
import polars as pl
import pyarrow as pa
import scipy.sparse as sp

N_GENES = 40

# probe, n_cells, n_background, sparse format, near-empty
PROBES = [
    ("A-A01", 120, 1200, "csr", False),
    ("A-A02", 80, 900, "csc", False),
    ("B-H12", 0, 600, "csr", False),
    ("C-D07", 40, 500, "csr", True),
]


def write_zst_tsv(path: str, df: pl.DataFrame) -> None:
    with pa.CompressedOutputStream(path, "zstd") as out:
        out.write(df.write_csv(separator="\t").encode())


def write_h5ad(path: str, barcodes: np.ndarray, counts: np.ndarray, probe: str, fmt: str):
    X = sp.csr_matrix(counts) if fmt == "csr" else sp.csc_matrix(counts)
    adata = ad.AnnData(
        X=X.astype(np.float32),
        obs=pd.DataFrame(index=[f"{b}-{probe}" for b in barcodes]),
        var=pd.DataFrame(index=[f"gene{i}" for i in range(counts.shape[1])]),
    )
    adata.write_h5ad(path)
