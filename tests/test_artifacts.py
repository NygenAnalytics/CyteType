import tempfile
from collections.abc import Iterable
from typing import Any

import duckdb
import pytest
import h5py
import anndata
import numpy as np
import pandas as pd
import scipy.sparse as sp
from pathlib import Path

from cytetype.core.artifacts import (
    save_features_matrix,
    save_obs_duckdb,
    _is_integer_valued,
)


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
    # artifacts.py binds WRITE_MEM_BUDGET at import, so patch it there.
    monkeypatch.setattr("cytetype.core.artifacts.WRITE_MEM_BUDGET", 256)
    progress_labels: list[str] = []

    def _recording_tqdm(iterable: Iterable[int], desc: str = "", **kwargs: Any) -> Any:
        progress_labels.append(desc)
        return iterable

    monkeypatch.setattr(
        "cytetype.core.artifacts._try_import_tqdm", lambda: _recording_tqdm
    )

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

    # The column-group scatter ran in several groups, not one.
    assert "Writing normalized counts (group 2)" in progress_labels

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


def test_save_features_matrix_hints_backed_mode_only_for_large_in_memory_data(
    tmp_path: Path,
    mock_adata: anndata.AnnData,
    monkeypatch: pytest.MonkeyPatch,
    log_messages: list[tuple[str, str]],
) -> None:
    def _backed_mode_hints() -> list[str]:
        return [message for _, message in log_messages if "backed mode" in message]

    save_features_matrix(
        out_file=str(tmp_path / "small.h5"), mat=mock_adata.X, col_batch=10
    )
    assert _backed_mode_hints() == []

    monkeypatch.setattr(
        "cytetype.core.artifacts.BACKED_MODE_HINT_MIN_CELLS", mock_adata.n_obs
    )
    save_features_matrix(
        out_file=str(tmp_path / "large.h5"), mat=mock_adata.X, col_batch=10
    )
    assert len(_backed_mode_hints()) == 1


def test_save_features_matrix_warns_with_cause_when_raw_group_fails(
    tmp_path: Path,
    mock_adata: anndata.AnnData,
    monkeypatch: pytest.MonkeyPatch,
    log_messages: list[tuple[str, str]],
) -> None:
    def _fail(*args: Any, **kwargs: Any) -> None:
        raise RuntimeError("raw counts exploded")

    monkeypatch.setattr("cytetype.core.artifacts._write_raw_group", _fail)
    out_path = tmp_path / "vars.h5"

    save_features_matrix(
        out_file=str(out_path),
        mat=mock_adata.X,
        raw_mat=mock_adata.X,
        col_batch=10,
    )

    assert (
        "WARNING",
        f"Failed to write raw counts group to {out_path}, skipping: raw counts exploded",
    ) in log_messages
    with h5py.File(out_path, "r") as f:
        assert "raw" not in f


def _obs_frame() -> pd.DataFrame:
    return pd.DataFrame(
        {"cluster": pd.Categorical(["a", "b", "a"]), "n_genes": [10, 20, 30]},
        index=["c1", "c2", "c3"],
    )


def test_save_obs_duckdb_spills_to_a_private_temporary_directory(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    system_tmp = tmp_path / "system-tmp"
    system_tmp.mkdir()
    monkeypatch.setattr(tempfile, "tempdir", str(system_tmp))
    spill_dirs: list[Path] = []
    connect = duckdb.connect

    def _connect(database: str, *args: Any, **kwargs: Any) -> Any:
        spill_dir = Path(kwargs["config"]["temp_directory"])
        assert spill_dir.is_dir()
        spill_dirs.append(spill_dir)
        return connect(database, *args, **kwargs)

    monkeypatch.setattr("cytetype.core.artifacts.duckdb.connect", _connect)
    out_path = tmp_path / "obs.duckdb"

    save_obs_duckdb(str(out_path), _obs_frame())

    assert len(spill_dirs) == 1
    assert spill_dirs[0].parent == system_tmp
    assert spill_dirs[0].name.startswith("cytetype-duckdb-")
    # Removed after writing.
    assert list(system_tmp.iterdir()) == []
    with connect(str(out_path), read_only=True) as con:
        assert con.execute("SELECT count(*) FROM obs").fetchone() == (3,)


def test_save_obs_duckdb_uses_and_keeps_an_explicit_temp_directory(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    spill_dir = tmp_path / "spill"
    seen: list[str] = []
    connect = duckdb.connect

    def _connect(database: str, *args: Any, **kwargs: Any) -> Any:
        seen.append(kwargs["config"]["temp_directory"])
        return connect(database, *args, **kwargs)

    monkeypatch.setattr("cytetype.core.artifacts.duckdb.connect", _connect)

    save_obs_duckdb(
        str(tmp_path / "obs.duckdb"), _obs_frame(), temp_directory=str(spill_dir)
    )

    assert seen == [str(spill_dir)]


def test_save_obs_duckdb_falls_back_when_temp_directory_cannot_be_created(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    log_messages: list[tuple[str, str]],
) -> None:
    def _no_temp_directory(*args: Any, **kwargs: Any) -> Any:
        raise OSError(28, "No space left on device")

    configs: list[dict[str, Any]] = []
    connect = duckdb.connect

    def _connect(database: str, *args: Any, **kwargs: Any) -> Any:
        configs.append(dict(kwargs["config"]))
        return connect(database, *args, **kwargs)

    monkeypatch.setattr(tempfile, "TemporaryDirectory", _no_temp_directory)
    monkeypatch.setattr("cytetype.core.artifacts.duckdb.connect", _connect)
    obs_df = _obs_frame()
    out_path = tmp_path / "obs.duckdb"

    save_obs_duckdb(
        str(out_path),
        obs_df,
        obsm_coordinates=np.arange(6, dtype=float).reshape(3, 2),
        coordinates_key="X_umap",
    )

    # DuckDB's default spill location is used, and obs_df is left as it was.
    assert "temp_directory" not in configs[0]
    assert list(obs_df.columns) == ["cluster", "n_genes"]
    assert any(
        level == "WARNING" and "No space left on device" in message
        for level, message in log_messages
    )
    with connect(str(out_path), read_only=True) as con:
        assert con.execute(
            'SELECT "__vis_coordinates_X_umap_2" FROM obs'
        ).fetchall() == [(1.0,), (3.0,), (5.0,)]


def test_save_obs_duckdb_removes_added_columns_when_writing_fails(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    def _connect(*args: Any, **kwargs: Any) -> Any:
        raise duckdb.IOException("disk full")

    monkeypatch.setattr("cytetype.core.artifacts.duckdb.connect", _connect)
    obs_df = _obs_frame()

    with pytest.raises(duckdb.IOException):
        save_obs_duckdb(
            str(tmp_path / "obs.duckdb"),
            obs_df,
            obsm_coordinates=np.zeros((3, 2)),
            coordinates_key="X_umap",
        )

    assert list(obs_df.columns) == ["cluster", "n_genes"]
