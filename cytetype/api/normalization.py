"""Normalization of per-cluster annotation results.

- Jobs since 2026-01-25 store the legacy fields (``justification``,
  ``supportingMarkers``, ...) empty and keep the content in
  ``latest.annotation.fullOutput.cellType``. Blank legacy fields are filled
  from there; non-blank values are never overwritten. The server's read-time
  fill in ``GET /results`` applies the same legacy-field rules, so both must
  change together.
- ``reviewerConfidence``, ``isHeterogeneous`` and ``heterogeneityType`` are
  added from ``latest.review`` (SDK-only keys).

Stored results may be partial or malformed, so every access is type-checked.
"""

from collections.abc import Mapping
from typing import Any

REVIEWER_CONFIDENCE_LEVELS: tuple[str, ...] = ("High", "Moderate", "Low")


def _as_dict(value: object) -> dict[str, Any]:
    return value if isinstance(value, dict) else {}


def _is_blank(value: object) -> bool:
    if value is None:
        return True
    if isinstance(value, str):
        return not value.strip()
    return isinstance(value, list) and not value


def _gene_list(value: object) -> list[str] | None:
    """Return the non-blank string items of ``value``, or None if it is not a list."""
    if not isinstance(value, list):
        return None
    return [item for item in value if isinstance(item, str) and item.strip()]


def _gene_text(reasoning: object, genes: object) -> str | None:
    if isinstance(reasoning, str) and reasoning.strip():
        return reasoning
    gene_list = _gene_list(genes)
    return None if gene_list is None else ", ".join(gene_list)


def _legacy_field_sources(
    cell_type: Mapping[str, Any],
) -> dict[str, str | list[str] | None]:
    """Map each legacy field to its value from ``fullOutput.cellType`` (None if unavailable)."""
    conclusion = cell_type.get("conclusion")
    return {
        "justification": conclusion if isinstance(conclusion, str) else None,
        "supportingMarkers": _gene_list(cell_type.get("keySupportingGenes")),
        "conflictingMarkers": _gene_list(cell_type.get("unexpectedGenes")),
        "missingExpression": _gene_text(
            cell_type.get("missingGenesReasoning"), cell_type.get("missingGenes")
        ),
        "unexpectedExpression": _gene_text(
            cell_type.get("unexpectedGenesReasoning"), cell_type.get("unexpectedGenes")
        ),
    }


def _review_fields(latest: Mapping[str, Any]) -> dict[str, Any]:
    review = _as_dict(latest.get("review"))
    confidence = review.get("confidence")
    is_heterogeneous = review.get("isHeterogeneous")
    heterogeneity_type = _as_dict(
        _as_dict(review.get("fullOutput")).get("heterogeneityAssessment")
    ).get("heterogeneityType")
    return {
        "reviewerConfidence": (
            confidence
            if isinstance(confidence, str) and confidence in REVIEWER_CONFIDENCE_LEVELS
            else None
        ),
        "isHeterogeneous": (
            is_heterogeneous if isinstance(is_heterogeneous, bool) else None
        ),
        "heterogeneityType": (
            heterogeneity_type
            if isinstance(heterogeneity_type, str) and heterogeneity_type.strip()
            else None
        ),
    }


def fill_cluster_annotation(
    entry: Mapping[str, Any],
    cluster_data: object,
) -> dict[str, Any]:
    """Return a normalized copy of one ``annotations`` entry (fill-only).

    Blank legacy fields are filled from ``latest.annotation.fullOutput.cellType``
    of the cluster's raw result; non-blank values are kept. The SDK-only keys are
    added when missing or None. Neither argument is modified.

    Args:
        entry: The cluster's entry in the ``annotations`` list.
        cluster_data: The cluster's raw result (``{"latest": ..., "history": ...}``).

    Returns:
        A new dict with the normalized entry.
    """
    normalized = dict(entry)
    latest = _as_dict(_as_dict(cluster_data).get("latest"))
    annotation = _as_dict(latest.get("annotation"))
    cell_type = _as_dict(annotation.get("fullOutput")).get("cellType")
    # A structured cellType marks gen-2 and V2 results; legacy results have none.
    if isinstance(cell_type, dict):
        for field, value in _legacy_field_sources(cell_type).items():
            if value is not None and _is_blank(normalized.get(field)):
                normalized[field] = value
    for key, value in _review_fields(latest).items():
        if normalized.get(key) is None:
            normalized[key] = value
    return normalized


def normalize_cluster_annotation(
    cluster_id: str,
    cluster_data: object,
) -> dict[str, Any] | None:
    """Build the normalized ``annotations`` entry for one raw cluster result.

    Args:
        cluster_id: The cluster's key in the API response's ``annotations``.
        cluster_data: The cluster's raw result (``{"latest": ..., "history": ...}``).

    Returns:
        The normalized entry, or None if the cluster has no annotation.
    """
    annotation = _as_dict(_as_dict(cluster_data).get("latest")).get("annotation")
    if not isinstance(annotation, dict):
        return None
    entry = {
        "clusterId": annotation.get("clusterId", cluster_id),
        "annotation": annotation.get("annotation", "Unknown"),
        "ontologyTerm": annotation.get("cellOntologyTermName", "Unknown"),
        "ontologyTermID": annotation.get("cellOntologyTerm", "Unknown"),
        "granularAnnotation": annotation.get("granularAnnotation", ""),
        "cellState": annotation.get("cellState", ""),
        "justification": annotation.get("justification", ""),
        "supportingMarkers": annotation.get("supportingMarkers", []),
        "conflictingMarkers": annotation.get("conflictingMarkers", []),
        "missingExpression": annotation.get("missingExpression", ""),
        "unexpectedExpression": annotation.get("unexpectedExpression", ""),
    }
    return fill_cluster_annotation(entry, cluster_data)


def renormalize_results(result: Mapping[str, Any]) -> dict[str, Any]:
    """Re-apply the fill-only normalization to stored results.

    Results saved by older versions gain the filled legacy fields and the
    SDK-only keys from their ``raw_annotations``. Without ``raw_annotations``
    (files from SDK 0.7.0 and older) the SDK-only keys are added as None, so
    every entry has the same keys.

    Args:
        result: Stored results, possibly the object held in ``adata.uns``.

    Returns:
        A new dict; ``result`` and its nested objects are not modified.
    """
    normalized = dict(result)
    annotations = result.get("annotations")
    if not isinstance(annotations, list):
        return normalized
    raw_annotations = result.get("raw_annotations")
    if not isinstance(raw_annotations, dict):
        raw_annotations = {}

    # fetch_job_results made one entry per raw cluster with an annotation dict,
    # in raw order, so match entries by consuming those clusters in order. This
    # keeps entries with a duplicated clusterId paired with their own cluster.
    pending: dict[str, list[Any]] = {}
    for cluster_id, cluster_data in raw_annotations.items():
        annotation = _as_dict(_as_dict(cluster_data).get("latest")).get("annotation")
        if not isinstance(annotation, dict):
            continue
        entry_id = annotation.get("clusterId", cluster_id)
        if isinstance(entry_id, str):
            pending.setdefault(entry_id, []).append(cluster_data)

    renormalized: list[Any] = []
    for entry in annotations:
        if not isinstance(entry, dict):
            renormalized.append(entry)
            continue
        cluster_id = entry.get("clusterId")
        cluster_data = None
        if isinstance(cluster_id, str):
            matches = pending.get(cluster_id)
            cluster_data = (
                matches.pop(0) if matches else raw_annotations.get(cluster_id)
            )
        renormalized.append(fill_cluster_annotation(entry, cluster_data))
    normalized["annotations"] = renormalized
    return normalized
