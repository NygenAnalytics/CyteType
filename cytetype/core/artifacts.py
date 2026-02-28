import re
from typing import Any, cast
from collections.abc import Sequence

import duckdb
import h5py
import hdf5plugin
import numpy as np
import pandas as pd
from pandas.api import types as ptypes
import scipy.sparse as sp
from anndata.abc import CSCDataset, CSRDataset

from ..config import logger


def _safe_column_dataset_name(
    source_name: str,
    column_index: int,
    group: h5py.Group,
) -> str:
    base = re.sub(r"[^A-Za-z0-9_.-]", "_", source_name).strip("_")
    if not base:
        base = f"column_{column_index}"
    candidate = base
    suffix = 1
    while candidate in group:
        candidate = f"{base}_{suffix}"
        suffix += 1
    return candidate


def _as_string_values(values: pd.Series | pd.Index | Sequence[Any]) -> np.ndarray:
    series = pd.Series(values, copy=False)
    return cast(np.ndarray, series.astype("string").fillna("").to_numpy(dtype=object))


def _write_var_metadata(
    out_file_group: h5py.File,
    n_cols: int,
    var_df: pd.DataFrame,
    var_names: pd.Index | Sequence[Any] | None,
) -> None:
    if len(var_df) != n_cols:
        raise ValueError(
            f"`adata.var` row count ({len(var_df)}) does not match matrix columns ({n_cols})."
        )

    names_source: pd.Index | Sequence[Any] = (
        var_names if var_names is not None else var_df.index
    )
    if len(names_source) != n_cols:
        raise ValueError(
            f"`var_names` length ({len(names_source)}) does not match matrix columns ({n_cols})."
        )

    info_group = out_file_group.create_group("info")
    var_group = info_group.create_group("var")
    text_dtype = h5py.string_dtype(encoding="utf-8")

    var_group.create_dataset(
        "var_names",
        data=_as_string_values(names_source),
        dtype=text_dtype,
    )
    var_group.create_dataset(
        "index",
        data=_as_string_values(var_df.index),
        dtype=text_dtype,
    )

    columns_group = var_group.create_group("columns")
    for i, col_name in enumerate(var_df.columns):
        source_name = str(col_name)
        dataset_name = _safe_column_dataset_name(source_name, i, columns_group)
        series = var_df[col_name]

        if isinstance(series.dtype, pd.CategoricalDtype):
            dataset = columns_group.create_dataset(
                dataset_name,
                data=_as_string_values(series.astype("string")),
                dtype=text_dtype,
            )
        elif ptypes.is_bool_dtype(series.dtype):
            if series.isna().any():
                bool_with_missing = series.astype("Int8").to_numpy(
                    dtype=np.int8,
                    na_value=-1,
                )
                dataset = columns_group.create_dataset(
                    dataset_name,
                    data=bool_with_missing,
                    dtype=np.int8,
                )
                dataset.attrs["missing_sentinel"] = -1
            else:
                dataset = columns_group.create_dataset(
                    dataset_name,
                    data=series.to_numpy(dtype=np.bool_),
                    dtype=np.bool_,
                )
        elif ptypes.is_numeric_dtype(series.dtype):
            numeric_data = pd.to_numeric(series, errors="coerce").to_numpy()
            dataset = columns_group.create_dataset(dataset_name, data=numeric_data)
        else:
            dataset = columns_group.create_dataset(
                dataset_name,
                data=_as_string_values(series),
                dtype=text_dtype,
            )

        dataset.attrs["source_name"] = source_name
        dataset.attrs["source_dtype"] = str(series.dtype)


def _is_integer_valued(mat: Any, sample_n_rows: int = 200) -> bool:
    if hasattr(mat, "dtype") and np.issubdtype(mat.dtype, np.integer):
        return True

    n_rows = mat.shape[0]
    row_end = min(sample_n_rows, n_rows)
    chunk = mat[:row_end]

    if sp.issparse(chunk):
        sample = chunk.data
    elif hasattr(chunk, "toarray"):
        sample = chunk.toarray().ravel()
    else:
        sample = np.asarray(chunk).ravel()

    if sample.size == 0:
        return True

    sample = sample.astype(np.float64, copy=False)
    return bool(np.all(np.isfinite(sample)) and np.all(sample == np.floor(sample)))


def _write_raw_group(
    f: h5py.File,
    mat: Any,
    n_obs: int,
    n_vars: int,
    col_indices: "np.ndarray | None",
    cell_batch: int,
    min_chunk_size: int,
) -> None:
    group = f.create_group("raw")
    group.attrs["n_obs"] = n_obs
    group.attrs["n_vars"] = n_vars

    chunk_size = max(1, min(n_vars * 10, min_chunk_size))

    d_indices = group.create_dataset(
        "indices",
        shape=(0,),
        maxshape=(n_obs * n_vars,),
        chunks=(chunk_size,),
        dtype=np.int32,
        compression=hdf5plugin.LZ4(),
    )
    d_data = group.create_dataset(
        "data",
        shape=(0,),
        maxshape=(n_obs * n_vars,),
        chunks=(chunk_size,),
        dtype=np.int32,
        compression=hdf5plugin.LZ4(),
    )

    indptr: list[int] = [0]
    for start in range(0, n_obs, cell_batch):
        end = min(start + cell_batch, n_obs)
        raw_chunk = mat[start:end]

        if col_indices is not None:
            raw_chunk = raw_chunk[:, col_indices]

        chunk = (
            raw_chunk.tocsr() if sp.issparse(raw_chunk) else sp.csr_matrix(raw_chunk)
        )

        old_size = d_indices.shape[0]
        chunk_col_indices = chunk.indices.astype(np.int32, copy=False)
        chunk_data = chunk.data.astype(np.int32, copy=False)

        d_indices.resize(old_size + len(chunk_col_indices), axis=0)
        d_indices[old_size : old_size + len(chunk_col_indices)] = chunk_col_indices

        d_data.resize(old_size + len(chunk_data), axis=0)
        d_data[old_size : old_size + len(chunk_data)] = chunk_data

        indptr.extend((chunk.indptr[1:] + indptr[-1]).tolist())

    group.create_dataset("indptr", data=np.asarray(indptr, dtype=np.int64))


def save_features_matrix(
    out_file: str,
    mat: Any,
    var_df: pd.DataFrame | None = None,
    var_names: pd.Index | Sequence[Any] | None = None,
    raw_mat: Any | None = None,
    raw_col_indices: "np.ndarray | None" = None,
    raw_cell_batch: int = 2000,
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

        if var_df is not None:
            _write_var_metadata(
                out_file_group=f,
                n_cols=n_cols,
                var_df=var_df,
                var_names=var_names,
            )

        if raw_mat is not None:
            try:
                _write_raw_group(
                    f,
                    raw_mat,
                    n_rows,
                    n_cols,
                    raw_col_indices,
                    raw_cell_batch,
                    min_chunk_size,
                )
            except Exception:
                logger.warning(
                    "Failed to write raw counts group to %s, skipping.", out_file
                )
                if "raw" in f:
                    del f["raw"]


def save_obs_duckdb(
    out_file: str,
    obs_df: pd.DataFrame,
    obsm_coordinates: np.ndarray | None = None,
    coordinates_key: str | None = None,
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

    if obsm_coordinates is not None and coordinates_key is not None:
        obs_df = obs_df.copy()
        obs_df[f"__vis_coordinates_{coordinates_key}_1"] = obsm_coordinates[:, 0]
        obs_df[f"__vis_coordinates_{coordinates_key}_2"] = obsm_coordinates[:, 1]

    dd_config: dict[str, Any] = {
        "threads": threads,
        "memory_limit": memory_limit,
        "temp_directory": temp_directory,
    }
    with duckdb.connect(out_file, config=dd_config) as con:
        con.register("obs_df", obs_df)
        con.execute(f"CREATE OR REPLACE TABLE {table_name} AS SELECT * FROM obs_df")
