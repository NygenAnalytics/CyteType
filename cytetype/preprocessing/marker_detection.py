from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Literal

import anndata
import numpy as np
import pandas as pd
from natsort import natsorted
from scipy.stats import ttest_ind_from_stats

from ..config import logger


@dataclass
class GroupStats:
    n: np.ndarray
    nnz: np.ndarray | None
    sum_: np.ndarray | None
    sum_sq: np.ndarray | None


def _accumulate_group_stats(
    X: Any,
    cell_group_indices: np.ndarray,
    n_groups: int,
    n_genes: int,
    cell_batch_size: int = 5000,
    compute_nnz: bool = False,
    compute_moments: bool = False,
    progress_desc: str = "streaming",
) -> GroupStats:
    """Single-pass row-batched accumulation of per-group statistics.

    Streams row chunks from *X* (works with dense, sparse, and backed
    matrices) and accumulates the requested statistics per group.
    """
    n_cells = X.shape[0]

    n_ = np.zeros(n_groups, dtype=np.int64)
    nnz_ = np.zeros((n_groups, n_genes), dtype=np.int64) if compute_nnz else None
    sum_ = np.zeros((n_groups, n_genes), dtype=np.float64) if compute_moments else None
    sum_sq_ = (
        np.zeros((n_groups, n_genes), dtype=np.float64) if compute_moments else None
    )

    chunk_starts = range(0, n_cells, cell_batch_size)
    try:
        import warnings

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            from tqdm.auto import tqdm

        chunk_iter = tqdm(chunk_starts, desc=progress_desc, unit="chunk")
    except ImportError:
        chunk_iter = chunk_starts

    for start in chunk_iter:
        end = min(start + cell_batch_size, n_cells)
        chunk = X[start:end]
        if hasattr(chunk, "toarray"):
            chunk = chunk.toarray()
        chunk = np.asarray(chunk, dtype=np.float64)
        chunk_labels = cell_group_indices[start:end]

        for g_idx in range(n_groups):
            mask = chunk_labels == g_idx
            if not mask.any():
                continue
            g_data = chunk[mask]
            n_[g_idx] += mask.sum()
            if sum_ is not None:
                sum_[g_idx] += g_data.sum(axis=0)
            if sum_sq_ is not None:
                sum_sq_[g_idx] += (g_data**2).sum(axis=0)
            if nnz_ is not None:
                nnz_[g_idx] += (g_data != 0).sum(axis=0)

    return GroupStats(n=n_, nnz=nnz_, sum_=sum_, sum_sq=sum_sq_)


def _benjamini_hochberg(pvals: np.ndarray) -> np.ndarray:
    n = len(pvals)
    order = np.argsort(pvals)
    pvals_sorted = pvals[order]
    adjusted = pvals_sorted * n / np.arange(1, n + 1)
    adjusted = np.minimum(adjusted, 1.0)
    # enforce monotonicity from the right (largest p-value first)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    result = np.empty_like(pvals)
    result[order] = adjusted
    return result


def rank_genes_groups_backed(
    adata: anndata.AnnData,
    groupby: str,
    *,
    use_raw: bool = False,
    layer: str | None = None,
    n_genes: int | None = None,
    key_added: str = "rank_genes_groups",
    cell_batch_size: int = 5000,
    corr_method: Literal["benjamini-hochberg", "bonferroni"] = "benjamini-hochberg",
    pts: bool = False,
    copy: bool = False,
) -> anndata.AnnData | None:
    """Memory-efficient rank_genes_groups for backed/on-disk AnnData objects.

    Drop-in replacement for ``sc.tl.rank_genes_groups`` that works on backed
    ``_CSRDataset`` matrices by streaming cell chunks instead of loading the
    full matrix. Uses Welch's t-test with one-vs-rest comparison.

    Writes scanpy-compatible output to ``adata.uns[key_added]``.
    """
    adata = adata.copy() if copy else adata

    # --- resolve data source ---
    if layer is not None:
        if use_raw:
            raise ValueError("Cannot specify `layer` and have `use_raw=True`.")
        X = adata.layers[layer]
        var_names = adata.var_names
    elif use_raw:
        if adata.raw is None:
            raise ValueError("Received `use_raw=True`, but `adata.raw` is empty.")
        X = adata.raw.X
        var_names = adata.raw.var_names
    else:
        X = adata.X
        var_names = adata.var_names

    n_cells, n_genes_total = X.shape

    # --- resolve groups ---
    groupby_series = adata.obs[groupby]
    if not hasattr(groupby_series, "cat"):
        groupby_series = groupby_series.astype("category")
    groups_order = np.array(
        natsorted(groupby_series.cat.categories.tolist()), dtype=object
    )
    group_labels = groupby_series.values

    # reject singlet groups
    group_counts = groupby_series.value_counts()
    singlets = group_counts[group_counts < 2].index.tolist()
    if singlets:
        raise ValueError(
            f"Could not calculate statistics for groups {', '.join(str(s) for s in singlets)} "
            "since they only contain one sample."
        )

    n_groups = len(groups_order)
    group_to_idx = {g: i for i, g in enumerate(groups_order)}
    cell_group_indices = np.array([group_to_idx[g] for g in group_labels])

    # --- log1p base handling (matches scanpy) ---
    log1p_base = adata.uns.get("log1p", {}).get("base")

    def expm1_func(x: np.ndarray) -> np.ndarray:
        if log1p_base is not None:
            result: np.ndarray = np.expm1(x * np.log(log1p_base))
            return result
        out: np.ndarray = np.expm1(x)
        return out

    # --- accumulate sufficient statistics in one pass ---
    stats = _accumulate_group_stats(
        X,
        cell_group_indices,
        n_groups,
        n_genes_total,
        cell_batch_size=cell_batch_size,
        compute_nnz=pts,
        compute_moments=True,
        progress_desc="rank_genes_groups_backed",
    )
    n_ = stats.n
    nnz_ = stats.nnz
    sum_ = stats.sum_
    sum_sq_ = stats.sum_sq
    assert sum_ is not None and sum_sq_ is not None
    total_sum = sum_.sum(axis=0)
    total_sum_sq = sum_sq_.sum(axis=0)

    # --- compute per-group statistics and t-test ---
    logger.info("Computing t-test statistics for {} groups...", n_groups)

    n_out = n_genes_total if n_genes is None else min(n_genes, n_genes_total)
    result_names = np.empty((n_out, n_groups), dtype=object)
    result_scores = np.empty((n_out, n_groups), dtype=np.float32)
    result_logfc = np.empty((n_out, n_groups), dtype=np.float32)
    result_pvals = np.empty((n_out, n_groups), dtype=np.float64)
    result_pvals_adj = np.empty((n_out, n_groups), dtype=np.float64)

    if pts and nnz_ is not None:
        pts_arr = np.zeros((n_groups, n_genes_total), dtype=np.float64)
        pts_rest_arr = np.zeros((n_groups, n_genes_total), dtype=np.float64)
        total_nnz = nnz_.sum(axis=0)

    for g_idx in range(n_groups):
        ng = n_[g_idx]
        nr = n_cells - ng

        mean_g = sum_[g_idx] / ng
        var_g = (sum_sq_[g_idx] - sum_[g_idx] ** 2 / ng) / (ng - 1)

        sum_rest = total_sum - sum_[g_idx]
        sum_sq_rest = total_sum_sq - sum_sq_[g_idx]
        mean_r = sum_rest / nr
        var_r = (sum_sq_rest - sum_rest**2 / nr) / (nr - 1)

        with np.errstate(invalid="ignore"):
            scores, pvals = ttest_ind_from_stats(
                mean1=mean_g,
                std1=np.sqrt(var_g),
                nobs1=ng,
                mean2=mean_r,
                std2=np.sqrt(var_r),
                nobs2=nr,
                equal_var=False,
            )

        scores[np.isnan(scores)] = 0
        pvals[np.isnan(pvals)] = 1

        if corr_method == "benjamini-hochberg":
            pvals_adj = _benjamini_hochberg(pvals)
        else:
            pvals_adj = np.minimum(pvals * n_genes_total, 1.0)

        with np.errstate(divide="ignore", invalid="ignore"):
            foldchanges = (expm1_func(mean_g) + 1e-9) / (expm1_func(mean_r) + 1e-9)
            logfc = np.log2(foldchanges)

        top_indices = np.argsort(-scores)[:n_out]

        result_names[:, g_idx] = var_names[top_indices]
        result_scores[:, g_idx] = scores[top_indices].astype(np.float32)
        result_logfc[:, g_idx] = logfc[top_indices].astype(np.float32)
        result_pvals[:, g_idx] = pvals[top_indices]
        result_pvals_adj[:, g_idx] = pvals_adj[top_indices]

        if pts and nnz_ is not None:
            pts_arr[g_idx] = nnz_[g_idx] / ng
            rest_nnz = total_nnz - nnz_[g_idx]
            pts_rest_arr[g_idx] = rest_nnz / nr

    # --- build scanpy-compatible recarrays ---
    group_names = [str(g) for g in groups_order]

    def _to_recarray(data: np.ndarray, dtype_str: str) -> np.recarray:
        df = pd.DataFrame(data, columns=group_names)
        rec: np.recarray = df.to_records(index=False, column_dtypes=dtype_str)
        return rec

    adata.uns[key_added] = {
        "params": {
            "groupby": groupby,
            "reference": "rest",
            "method": "t-test",
            "use_raw": use_raw,
            "layer": layer,
            "corr_method": corr_method,
        },
        "names": _to_recarray(result_names, "O"),
        "scores": _to_recarray(result_scores, "float32"),
        "logfoldchanges": _to_recarray(result_logfc, "float32"),
        "pvals": _to_recarray(result_pvals, "float64"),
        "pvals_adj": _to_recarray(result_pvals_adj, "float64"),
    }

    if pts and nnz_ is not None:
        adata.uns[key_added]["pts"] = pd.DataFrame(
            pts_arr.T, index=var_names, columns=group_names
        )
        adata.uns[key_added]["pts_rest"] = pd.DataFrame(
            pts_rest_arr.T, index=var_names, columns=group_names
        )

    logger.info(
        "rank_genes_groups_backed complete — {} genes per group written to adata.uns['{}']",
        n_out,
        key_added,
    )

    if copy:
        return adata
    return None
