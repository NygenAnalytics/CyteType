import pytest
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
        gene_symbols_column="gene_symbols",
        col_batch=10,
    )

    with h5py.File(out_path, "r") as f:
        assert "vars" in f
        assert "info" in f
        assert "var" in f["info"]
        assert "columns" in f["info/var"]
        assert len(f["info/var/var_names"]) == mock_adata.n_vars
        assert len(f["info/var/index"]) == mock_adata.n_vars
        assert f["info/var"].attrs["gene_symbols_column"] == "gene_symbols"

        columns_group = f["info/var/columns"]
        assert len(columns_group.keys()) == mock_adata.var.shape[1]
        for dataset in columns_group.values():
            assert "source_name" in dataset.attrs
            assert "source_dtype" in dataset.attrs


def test_save_features_matrix_omits_gene_symbols_attr_when_not_provided(
    tmp_path: Path,
    mock_adata: anndata.AnnData,
) -> None:
    out_path = tmp_path / "vars.h5"

    save_features_matrix(
        out_file=str(out_path),
        mat=mock_adata.X,
        var_df=mock_adata.var,
        var_names=mock_adata.var_names,
        gene_symbols_column=None,
        col_batch=10,
    )

    with h5py.File(out_path, "r") as f:
        assert "gene_symbols_column" not in f["info/var"].attrs


def test_save_features_matrix_omits_gene_symbols_attr_when_omitted(
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
        assert "gene_symbols_column" not in f["info/var"].attrs


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


def test_save_features_matrix_backed_csr(tmp_path: Path) -> None:
    n_obs, n_vars = 200, 80
    rng = np.random.default_rng(42)
    mat = sp.random(n_obs, n_vars, density=0.3, format="csr", random_state=rng)
    mat.data = rng.standard_normal(mat.nnz).astype(np.float32)
    reference_csc = mat.tocsc()

    h5ad_path = tmp_path / "backed.h5ad"
    adata = anndata.AnnData(X=mat)
    adata.write_h5ad(h5ad_path)
    del adata

    backed = anndata.read_h5ad(h5ad_path, backed="r")
    out_path = tmp_path / "vars.h5"
    save_features_matrix(out_file=str(out_path), mat=backed.X)
    backed.file.close()

    with h5py.File(out_path, "r") as f:
        grp = f["vars"]
        assert grp.attrs["n_obs"] == n_obs
        assert grp.attrs["n_vars"] == n_vars

        h5_indices = grp["indices"][:]
        h5_data = grp["data"][:]
        h5_indptr = grp["indptr"][:]

        assert h5_indices.dtype == np.int32
        assert h5_data.dtype == np.float32
        assert len(h5_indptr) == n_vars + 1
        assert h5_indptr[0] == 0
        assert h5_indptr[-1] == reference_csc.nnz

        np.testing.assert_array_equal(h5_indptr, reference_csc.indptr)
        np.testing.assert_array_equal(h5_indices, reference_csc.indices)
        np.testing.assert_allclose(h5_data, reference_csc.data, rtol=1e-6)


def test_save_features_matrix_backed_csr_multiple_groups(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    import cytetype.config

    monkeypatch.setattr(cytetype.config, "WRITE_MEM_BUDGET", 256)

    n_obs, n_vars = 50, 30
    rng = np.random.default_rng(7)
    mat = sp.random(n_obs, n_vars, density=0.4, format="csr", random_state=rng)
    mat.data = rng.standard_normal(mat.nnz).astype(np.float32)
    reference_csc = mat.tocsc()

    h5ad_path = tmp_path / "backed.h5ad"
    adata = anndata.AnnData(X=mat)
    adata.write_h5ad(h5ad_path)
    del adata

    backed = anndata.read_h5ad(h5ad_path, backed="r")
    out_path = tmp_path / "vars.h5"
    save_features_matrix(out_file=str(out_path), mat=backed.X)
    backed.file.close()

    with h5py.File(out_path, "r") as f:
        grp = f["vars"]
        h5_indptr = grp["indptr"][:]
        h5_indices = grp["indices"][:]
        h5_data = grp["data"][:]

        np.testing.assert_array_equal(h5_indptr, reference_csc.indptr)
        np.testing.assert_array_equal(h5_indices, reference_csc.indices)
        np.testing.assert_allclose(h5_data, reference_csc.data, rtol=1e-6)


def test_is_integer_valued_with_true_integers() -> None:
    mat = sp.csr_matrix(np.array([[1, 0, 3], [0, 2, 0]], dtype=np.int32))
    assert _is_integer_valued(mat) is True


def test_is_integer_valued_with_float_integers() -> None:
    mat = sp.csr_matrix(np.array([[1.0, 0.0, 3.0], [0.0, 2.0, 0.0]], dtype=np.float32))
    assert _is_integer_valued(mat) is True


def test_is_integer_valued_with_floats() -> None:
    mat = sp.csr_matrix(np.array([[1.5, 0.0, 3.2], [0.0, 2.7, 0.0]], dtype=np.float32))
    assert _is_integer_valued(mat) is False
