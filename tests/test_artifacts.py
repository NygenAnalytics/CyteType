import h5py
import anndata
import numpy as np
import scipy.sparse as sp
from pathlib import Path

from cytetype.core.artifacts import save_features_matrix, _is_integer_valued


def test_save_features_matrix_writes_var_metadata(
    tmp_path: Path,
    mock_adata: anndata.AnnData,
) -> None:
    out_path = tmp_path / "vars.h5"

    save_features_matrix(
        out_file=str(out_path),
        mat=mock_adata.X,
        var_df=mock_adata.var,
        var_names=mock_adata.var_names,
        col_batch=10,
    )

    with h5py.File(out_path, "r") as f:
        assert "vars" in f
        assert "info" in f
        assert "var" in f["info"]
        assert "columns" in f["info/var"]
        assert len(f["info/var/var_names"]) == mock_adata.n_vars
        assert len(f["info/var/index"]) == mock_adata.n_vars

        columns_group = f["info/var/columns"]
        assert len(columns_group.keys()) == mock_adata.var.shape[1]
        for dataset in columns_group.values():
            assert "source_name" in dataset.attrs
            assert "source_dtype" in dataset.attrs


def test_save_features_matrix_writes_raw_group(
    tmp_path: Path,
    mock_adata: anndata.AnnData,
) -> None:
    n_obs, n_vars = mock_adata.n_obs, mock_adata.n_vars
    rng = np.random.default_rng(0)
    raw_mat = sp.random(n_obs, n_vars, density=0.1, format="csr", random_state=rng)
    raw_mat.data = rng.integers(1, 20, size=raw_mat.nnz).astype(np.int32)

    out_path = tmp_path / "vars.h5"
    save_features_matrix(
        out_file=str(out_path),
        mat=mock_adata.X,
        var_df=mock_adata.var,
        var_names=mock_adata.var_names,
        raw_mat=raw_mat,
        raw_cell_batch=30,
        col_batch=10,
    )

    with h5py.File(out_path, "r") as f:
        assert "raw" in f
        raw = f["raw"]
        assert raw.attrs["n_obs"] == n_obs
        assert raw.attrs["n_vars"] == n_vars
        assert raw["data"].dtype == np.int32
        assert raw["indices"].dtype == np.int32
        assert len(raw["indptr"]) == n_obs + 1
        assert raw["indptr"][0] == 0
        assert raw["indptr"][-1] == raw_mat.nnz


def test_save_features_matrix_raw_skipped_when_none(
    tmp_path: Path,
    mock_adata: anndata.AnnData,
) -> None:
    out_path = tmp_path / "vars.h5"
    save_features_matrix(
        out_file=str(out_path),
        mat=mock_adata.X,
        col_batch=10,
    )

    with h5py.File(out_path, "r") as f:
        assert "raw" not in f
        assert "vars" in f


def test_save_features_matrix_raw_with_float_integers(
    tmp_path: Path,
    mock_adata: anndata.AnnData,
) -> None:
    n_obs, n_vars = mock_adata.n_obs, mock_adata.n_vars
    rng = np.random.default_rng(1)
    int_vals = rng.integers(1, 10, size=(n_obs, n_vars))
    raw_mat = sp.csr_matrix(int_vals.astype(np.float32))

    out_path = tmp_path / "vars.h5"
    save_features_matrix(
        out_file=str(out_path),
        mat=mock_adata.X,
        raw_mat=raw_mat,
        raw_cell_batch=30,
        col_batch=10,
    )

    with h5py.File(out_path, "r") as f:
        assert "raw" in f
        assert f["raw"]["data"].dtype == np.int32


def test_is_integer_valued_with_true_integers() -> None:
    mat = sp.csr_matrix(np.array([[1, 0, 3], [0, 2, 0]], dtype=np.int32))
    assert _is_integer_valued(mat) is True


def test_is_integer_valued_with_float_integers() -> None:
    mat = sp.csr_matrix(np.array([[1.0, 0.0, 3.0], [0.0, 2.0, 0.0]], dtype=np.float32))
    assert _is_integer_valued(mat) is True


def test_is_integer_valued_with_floats() -> None:
    mat = sp.csr_matrix(np.array([[1.5, 0.0, 3.2], [0.0, 2.7, 0.0]], dtype=np.float32))
    assert _is_integer_valued(mat) is False
