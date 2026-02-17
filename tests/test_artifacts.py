import h5py
import anndata
from pathlib import Path

from cytetype.core.artifacts import save_features_matrix


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
