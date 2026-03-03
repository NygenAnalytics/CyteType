from __future__ import annotations

import anndata
import pandas as pd
from natsort import natsorted

from ..config import logger


def subsample_by_group(
    adata: anndata.AnnData,
    group_key: str,
    max_cells_per_group: int = 1000,
    random_state: int = 0,
) -> anndata.AnnData:
    """Subsample cells from an AnnData object, capping each group to a maximum count.

    Groups smaller than *max_cells_per_group* are kept intact.

    Parameters
    ----------
    adata
        The AnnData object to subsample.
    group_key
        Column in ``adata.obs`` that defines the groups (e.g. cluster labels).
    max_cells_per_group
        Maximum number of cells to retain per group. Groups with fewer cells
        are included in full.
    random_state
        Seed for reproducible sampling.

    Returns
    -------
    anndata.AnnData
        A new in-memory AnnData object containing at most *max_cells_per_group*
        cells per group.
    """
    if group_key not in adata.obs.columns:
        raise KeyError(
            f"Group key '{group_key}' not found in adata.obs. "
            f"Available columns: {list(adata.obs.columns)}"
        )

    is_backed = getattr(adata, "isbacked", False)
    groups = natsorted(adata.obs[group_key].unique())
    subsampled: list[anndata.AnnData] = []

    try:
        import warnings

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            from tqdm.auto import tqdm

        group_iter = tqdm(groups, desc="Subsampling groups", unit="group")
    except ImportError:
        group_iter = groups

    for group in group_iter:
        mask = adata.obs[group_key] == group
        n_cells = mask.sum()

        if n_cells > max_cells_per_group:
            keep = pd.Series(False, index=adata.obs.index)
            sampled = mask[mask].sample(n=max_cells_per_group, random_state=random_state)
            keep[sampled.index] = True
            subset = adata[keep]
        else:
            subset = adata[mask]

        subsampled.append(subset.to_memory() if is_backed else subset.copy())

    result = anndata.concat(subsampled)

    logger.info(f"Subsampling complete: {adata.n_obs} -> {result.n_obs} cells")

    return result
