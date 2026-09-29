"""Constants and writers for the synthetic cyto runs in ``conftest.py``."""

import polars as pl
import pyarrow as pa

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
