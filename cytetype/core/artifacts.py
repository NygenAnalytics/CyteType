import re
from typing import Any

import duckdb
import h5py
import hdf5plugin
import numpy as np
import pandas as pd
import scipy.sparse as sp
from anndata.abc import CSCDataset, CSRDataset

from ..config import logger


def save_features_matrix(
    out_file: str,
    mat: Any,
    min_chunk_size: int = 10_000_000,
    col_batch: int | None = None,
) -> None:
    """Persist normalized count matrix in the vars.h5 format expected by the API."""
    if not hasattr(mat, "shape") or not hasattr(mat, "__getitem__"):
        raise TypeError(
            "mat must be an array-like matrix with shape and slicing support"
        )

    n_rows, n_cols = mat.shape

    if not isinstance(mat, (CSRDataset, CSCDataset)):
        logger.warning(
            "For large datasets, use AnnData backed mode (e.g. sc.read_h5ad(..., backed='r')) "
            "so `adata.X` is a backed sparse dataset and avoids loading the full matrix in memory.",
        )

    if col_batch is None:
        col_batch = max(1, int(100_000_000 / max(n_rows, 1)))

    chunk_size = max(1, min(n_rows * 10, min_chunk_size))
    with h5py.File(out_file, "w") as f:
        group = f.create_group("vars")
        group.attrs["n_obs"] = n_rows
        group.attrs["n_vars"] = n_cols

        d_indices = group.create_dataset(
            "indices",
            shape=(0,),
            maxshape=(n_rows * n_cols,),
            chunks=(chunk_size,),
            dtype=np.int32,
            compression=hdf5plugin.LZ4(),
        )
        d_data = group.create_dataset(
            "data",
            shape=(0,),
            maxshape=(n_rows * n_cols,),
            chunks=(chunk_size,),
            dtype=np.float32,
            compression=hdf5plugin.LZ4(),
        )

        indptr: list[int] = [0]
        for start in range(0, n_cols, col_batch):
            end = min(start + col_batch, n_cols)
            raw_chunk = mat[:, start:end]
            chunk = (
                raw_chunk.tocsc()
                if sp.issparse(raw_chunk)
                else sp.csc_matrix(raw_chunk)
            )

            old_size = d_indices.shape[0]
            chunk_indices = chunk.indices.astype(np.int32, copy=False)
            chunk_data = chunk.data.astype(np.float32, copy=False)

            d_indices.resize(old_size + len(chunk_indices), axis=0)
            d_indices[old_size : old_size + len(chunk_indices)] = chunk_indices

            d_data.resize(old_size + len(chunk_data), axis=0)
            d_data[old_size : old_size + len(chunk_data)] = chunk_data

            indptr.extend((chunk.indptr[1:] + indptr[-1]).tolist())

        group.create_dataset("indptr", data=np.asarray(indptr, dtype=np.int64))


def save_obs_duckdb(
    out_file: str,
    obs_df: pd.DataFrame,
    table_name: str = "obs",
    threads: int = 4,
    memory_limit: str = "4GB",
    temp_directory: str = "/tmp/duckdb",
) -> None:
    """Persist adata.obs to a DuckDB table."""
    if not re.fullmatch(r"[A-Za-z_][A-Za-z0-9_]*", table_name):
        raise ValueError(
            "Invalid table_name. Use letters, numbers, and underscores only."
        )

    dd_config: dict[str, Any] = {
        "threads": threads,
        "memory_limit": memory_limit,
        "temp_directory": temp_directory,
    }
    with duckdb.connect(out_file, config=dd_config) as con:
        con.register("obs_df", obs_df)
        con.execute(f"CREATE OR REPLACE TABLE {table_name} AS SELECT * FROM obs_df")
