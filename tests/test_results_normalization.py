"""Normalization of /results annotations and how results are stored and read."""

import copy
import json
from pathlib import Path
from typing import Any

import anndata
import numpy as np
import pandas as pd
import pytest

from cytetype.api import client
from cytetype.api.normalization import (
    fill_cluster_annotation,
    normalize_cluster_annotation,
    renormalize_results,
)
from cytetype.api.transport import HTTPTransport
from cytetype.core.results import load_local_results, store_annotations
from cytetype.plotting import marker_dotplot

LEGACY_FIELDS = (
    "justification",
    "supportingMarkers",
    "conflictingMarkers",
    "missingExpression",
    "unexpectedExpression",
)
EXISTING_KEYS = (
    "clusterId",
    "annotation",
    "ontologyTerm",
    "ontologyTermID",
    "granularAnnotation",
    "cellState",
    *LEGACY_FIELDS,
)
SDK_KEYS = ("reviewerConfidence", "isHeterogeneous", "heterogeneityType")
EMPTY_LEGACY_FIELDS: dict[str, Any] = {
    "justification": "",
    "supportingMarkers": [],
    "conflictingMarkers": [],
    "missingExpression": "",
    "unexpectedExpression": "",
}
CONCLUSION = (
    "**CD8+ T cell**: `CD8A` and `CD8B` with cytotoxic `GZMK` support a\n"
    "cytotoxic T cell identity. "
)
EXPECTED_FILL: dict[str, Any] = {
    "justification": CONCLUSION,
    "supportingMarkers": ["CD8A", "CD8B", "GZMK"],
    "conflictingMarkers": ["MKI67"],
    "missingExpression": "`CCR7` and `SELL` mark naive cells.",
    "unexpectedExpression": "`MKI67` suggests a proliferating subset.",
}
LEGACY_EVIDENCE: dict[str, Any] = {
    "justification": "Strong GFAP and S100B expression.",
    "supportingMarkers": ["GFAP", "S100B", "AQP4"],
    "conflictingMarkers": [],
    "missingExpression": "SOX9 lower than expected",
    "unexpectedExpression": "",
}


def _gen2_cell_type(**overrides: Any) -> dict[str, Any]:
    """``fullOutput.cellType`` as stored between January and August 2026."""
    cell_type: dict[str, Any] = {
        "label": "Cd8+ T Cell",
        "keySupportingGenes": ["CD8A", "CD8B", "GZMK"],
        "keySupportingGenesReasoning": "`CD8A` and `CD8B` are lineage markers.",
        "missingGenes": ["CCR7", "SELL"],
        "missingGenesReasoning": "`CCR7` and `SELL` mark naive cells.",
        "unexpectedGenes": ["MKI67"],
        "unexpectedGenesReasoning": "`MKI67` suggests a proliferating subset.",
        "pathwaysSupportingEvidence": ["T cell receptor signaling"],
        "pathwaysSupportingEvidenceReasoning": "TCR signaling is enriched.",
        "candidateComparison": ["NK cell: lacks `NCAM1`."],
        "tissueContextFit": "Expected in blood.",
        "clusterContextFit": "Present in all samples.",
        "conclusion": CONCLUSION,
        "citationExplanations": {},
    }
    cell_type.update(overrides)
    return cell_type


def _v2_cell_type(**overrides: Any) -> dict[str, Any]:
    """``fullOutput.cellType`` as stored since August 2026."""
    cell_type = _gen2_cell_type(
        catalogReferenceLabel="CD8-positive, alpha-beta T cell",
        classification="catalogue",
        whyNotDeeper="",
        annotationLevels={
            "L1": {"label": "Immune cell"},
            "L2": {"label": "T cell"},
            "L3": {"label": "CD8+ T cell"},
        },
        candidateComparison={"NK cell": {"verdict": "rejected"}},
    )
    del cell_type["citationExplanations"]
    cell_type.update(overrides)
    return cell_type


def _review(
    confidence: Any,
    is_heterogeneous: Any,
    heterogeneity_type: Any = None,
) -> dict[str, Any]:
    return {
        "reviewComments": "",
        "confidence": confidence,
        "isHeterogeneous": is_heterogeneous,
        "feedback": "",
        "labelComparisonScore": None,
        "fullOutput": {
            "confidenceAssessment": {"label": confidence, "reasoning": ""},
            "heterogeneityAssessment": {
                "label": "Heterogeneous" if is_heterogeneous else "Homogeneous",
                "reasoning": "",
                "heterogeneityType": heterogeneity_type,
                "heterogeneityTypeReasoning": "",
            },
        },
    }


def _new_result(
    cluster_id: str,
    cell_type: Any,
    review: dict[str, Any] | None = None,
    label: str = "CD8+ T Cell",
    **legacy_fields: Any,
) -> dict[str, Any]:
    """A ``latest`` entry written by the orchestrator since January 2026."""
    return {
        "context": "Cluster context",
        "annotation": {
            "clusterId": cluster_id,
            "annotation": label,
            "granularAnnotation": "GZMK+ CD8+ T cell",
            "cellState": "Cytotoxic",
            "cellOntologyTerm": "CL:0000625",
            "cellOntologyTermName": "CD8-positive, alpha-beta T cell",
            **EMPTY_LEGACY_FIELDS,
            **legacy_fields,
            "fullOutput": {
                "clusterId": cluster_id,
                "cellType": cell_type,
                "cellStateWithReasoning": {"label": "Cytotoxic", "reasoning": ""},
            },
            "clusterContext": {},
        },
        "review": review if review is not None else _review("High", False),
        "developmentalSimilarity": [],
        "functionalSimilarity": [],
    }


def _legacy_result(cluster_id: str, **annotation_fields: Any) -> dict[str, Any]:
    """A ``latest`` entry stored before January 2026: evidence in the legacy fields."""
    return {
        "context": "Cluster context",
        "annotation": {
            "clusterId": cluster_id,
            "annotation": "Astrocyte",
            "granularAnnotation": "GFAP+ astrocyte",
            "cellState": "mature",
            "cellOntologyTerm": "CL:0000127",
            "cellOntologyTermName": "astrocyte",
            **LEGACY_EVIDENCE,
            "fullOutput": None,
            "clusterContext": None,
            **annotation_fields,
        },
        "review": {"confidence": "High"},
        "developmentalSimilarity": [],
        "functionalSimilarity": [],
    }


def _results_response() -> dict[str, Any]:
    """A ``GET /results`` response mixing every stored result generation."""
    return {
        "studyContext": "Human PBMC",
        "annotations": {
            "1": {"latest": _legacy_result("1"), "history": []},
            "2": {
                "latest": _new_result(
                    "2", _gen2_cell_type(), _review("Moderate", False)
                ),
                "history": [],
            },
            "3": {
                "latest": _new_result(
                    "3", _v2_cell_type(), _review("Low", True, "distinct_populations")
                ),
                "history": [_legacy_result("3")],
            },
            # A cluster without a result, and an older test-fixture shape.
            "4": {
                "latest": {"developmentalSimilarity": [], "functionalSimilarity": []},
                "history": [],
            },
            "5": {"latest": {"annotation": "B cell"}, "history": []},
        },
        "summary": {"tissue": "blood"},
        "semanticOrder": ["1", "2", "3"],
        "clusterCategories": [{"categoryName": "Lymphoid", "clusterIds": ["2", "3"]}],
        "displayLevel": "",
    }


def _fetch(monkeypatch: pytest.MonkeyPatch, response: dict[str, Any]) -> dict[str, Any]:
    requested: list[str] = []

    def _get(
        self: HTTPTransport, endpoint: str, timeout: int = 30
    ) -> tuple[int, dict[str, Any]]:
        requested.append(endpoint)
        return 200, response

    monkeypatch.setattr(HTTPTransport, "get", _get)
    result = client.fetch_job_results("https://api.test", "token", "job-1")
    assert requested == ["results/job-1"]
    return result


def _old_sdk_result(response: dict[str, Any]) -> dict[str, Any]:
    """What SDK versions before the normalization stored for ``response``."""
    annotations = []
    for cluster_id, cluster_data in response["annotations"].items():
        annotation = cluster_data["latest"].get("annotation")
        if isinstance(annotation, dict):
            annotations.append(
                {
                    "clusterId": annotation.get("clusterId", cluster_id),
                    "annotation": annotation.get("annotation", "Unknown"),
                    "ontologyTerm": annotation.get("cellOntologyTermName", "Unknown"),
                    "ontologyTermID": annotation.get("cellOntologyTerm", "Unknown"),
                    "granularAnnotation": annotation.get("granularAnnotation", ""),
                    "cellState": annotation.get("cellState", ""),
                    **{
                        field: annotation.get(field, default)
                        for field, default in EMPTY_LEGACY_FIELDS.items()
                    },
                }
            )
    return {
        "annotations": annotations,
        "summary": response["summary"],
        "clusterCategories": response["clusterCategories"],
        "studyContext": response["studyContext"],
        "raw_annotations": response["annotations"],
    }


def _legacy_values(entry: dict[str, Any]) -> dict[str, Any]:
    return {field: entry[field] for field in LEGACY_FIELDS}


def _sdk_values(entry: dict[str, Any]) -> dict[str, Any]:
    return {key: entry[key] for key in SDK_KEYS}


def _by_cluster(result: dict[str, Any]) -> dict[str, dict[str, Any]]:
    return {entry["clusterId"]: entry for entry in result["annotations"]}


def _normalize(latest: dict[str, Any]) -> dict[str, Any]:
    entry = normalize_cluster_annotation("1", {"latest": latest, "history": []})
    assert entry is not None
    return entry


def test_fetch_job_results_normalizes_every_result_generation(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    response = _results_response()
    before = copy.deepcopy(response)

    result = _fetch(monkeypatch, response)
    annotations = _by_cluster(result)

    # Clusters without an annotation dict are skipped, as before.
    assert list(annotations) == ["1", "2", "3"]
    assert _legacy_values(annotations["1"]) == LEGACY_EVIDENCE
    assert _legacy_values(annotations["2"]) == EXPECTED_FILL
    assert _legacy_values(annotations["3"]) == EXPECTED_FILL
    assert _sdk_values(annotations["1"]) == {
        "reviewerConfidence": "High",
        "isHeterogeneous": None,
        "heterogeneityType": None,
    }
    assert _sdk_values(annotations["2"]) == {
        "reviewerConfidence": "Moderate",
        "isHeterogeneous": False,
        "heterogeneityType": None,
    }
    assert _sdk_values(annotations["3"]) == {
        "reviewerConfidence": "Low",
        "isHeterogeneous": True,
        "heterogeneityType": "distinct_populations",
    }

    # Existing keys keep their names, order and values; the new keys come last.
    assert list(annotations["3"]) == [*EXISTING_KEYS, *SDK_KEYS]
    assert {key: annotations["3"][key] for key in EXISTING_KEYS[:6]} == {
        "clusterId": "3",
        "annotation": "CD8+ T Cell",
        "ontologyTerm": "CD8-positive, alpha-beta T cell",
        "ontologyTermID": "CL:0000625",
        "granularAnnotation": "GZMK+ CD8+ T cell",
        "cellState": "Cytotoxic",
    }

    # The raw API data is returned and left as received.
    assert result["raw_annotations"] == before["annotations"]
    assert response == before
    assert result["summary"] == before["summary"]
    assert result["clusterCategories"] == before["clusterCategories"]
    assert result["studyContext"] == before["studyContext"]


def test_non_blank_fields_are_never_overwritten() -> None:
    curated: dict[str, Any] = {
        "justification": "Curated justification",
        "supportingMarkers": ["CD3E"],
        "conflictingMarkers": ["CD19"],
        "missingExpression": "Curated missing note",
        "unexpectedExpression": "Curated unexpected note",
    }

    entry = _normalize(_new_result("1", _v2_cell_type(), **curated))

    assert _legacy_values(entry) == curated


def test_only_blank_fields_are_filled() -> None:
    entry = _normalize(
        _new_result(
            "1",
            _v2_cell_type(),
            justification="Curated justification",
            supportingMarkers=None,
            conflictingMarkers=["CD19"],
            missingExpression="   ",
            unexpectedExpression="\n",
        )
    )

    assert _legacy_values(entry) == {
        "justification": "Curated justification",
        "supportingMarkers": ["CD8A", "CD8B", "GZMK"],
        "conflictingMarkers": ["CD19"],
        "missingExpression": "`CCR7` and `SELL` mark naive cells.",
        "unexpectedExpression": "`MKI67` suggests a proliferating subset.",
    }


def test_missing_legacy_fields_are_filled_when_a_source_exists() -> None:
    latest = _new_result("1", _gen2_cell_type())
    for field in LEGACY_FIELDS:
        del latest["annotation"][field]

    assert _legacy_values(_normalize(latest)) == EXPECTED_FILL


def test_text_is_copied_as_stored() -> None:
    reasoning = "  `CCR7` is *absent*.\n"

    entry = _normalize(_new_result("1", _v2_cell_type(missingGenesReasoning=reasoning)))

    assert entry["justification"] == CONCLUSION
    assert entry["missingExpression"] == reasoning


@pytest.mark.parametrize("reasoning", ["", "  \n ", None, 7, ["CCR7"]])
def test_blank_reasoning_falls_back_to_the_joined_gene_list(reasoning: Any) -> None:
    entry = _normalize(
        _new_result(
            "1",
            _gen2_cell_type(
                missingGenesReasoning=reasoning, unexpectedGenesReasoning=reasoning
            ),
        )
    )

    assert entry["missingExpression"] == "CCR7, SELL"
    assert entry["unexpectedExpression"] == "MKI67"


def test_gene_lists_keep_only_non_blank_strings() -> None:
    noisy = [" CD8A", None, 5, "", "   ", {"gene": "CD8B"}, ["GZMK"], True, "GZMK"]

    entry = _normalize(
        _new_result(
            "1",
            _v2_cell_type(
                keySupportingGenes=noisy,
                unexpectedGenes=noisy,
                unexpectedGenesReasoning="",
                missingGenes=noisy,
                missingGenesReasoning="",
            ),
        )
    )

    assert entry["supportingMarkers"] == [" CD8A", "GZMK"]
    assert entry["conflictingMarkers"] == [" CD8A", "GZMK"]
    assert entry["missingExpression"] == " CD8A, GZMK"
    assert entry["unexpectedExpression"] == " CD8A, GZMK"


def test_empty_sources_fill_typed_blank_values_and_never_none() -> None:
    entry = _normalize(
        _new_result(
            "1",
            _v2_cell_type(
                conclusion="",
                keySupportingGenes=[None, ""],
                unexpectedGenes=[],
                unexpectedGenesReasoning="",
                missingGenes=[],
                missingGenesReasoning="",
            ),
            justification=None,
            supportingMarkers=None,
            conflictingMarkers=None,
            missingExpression=None,
            unexpectedExpression=None,
        )
    )

    assert _legacy_values(entry) == EMPTY_LEGACY_FIELDS


def test_missing_or_wrongly_typed_sources_leave_fields_unchanged() -> None:
    cell_type = {
        "conclusion": ["not", "a", "string"],
        "keySupportingGenes": "CD8A, CD8B",
        "unexpectedGenes": {"MKI67": 1},
        "missingGenesReasoning": "   ",
        "unexpectedGenesReasoning": None,
    }

    entry = _normalize(_new_result("1", cell_type, justification=None))

    assert _legacy_values(entry) == {**EMPTY_LEGACY_FIELDS, "justification": None}


@pytest.mark.parametrize(
    "full_output",
    [None, "not a dict", [], {}, {"cellType": "T cell"}, {"cellType": ["T cell"]}],
)
def test_results_without_a_structured_cell_type_keep_their_fields(
    full_output: Any,
) -> None:
    latest = _legacy_result("1", fullOutput=full_output, justification="")

    entry = _normalize(latest)

    assert _legacy_values(entry) == {**LEGACY_EVIDENCE, "justification": ""}


@pytest.mark.parametrize(
    ("review", "expected"),
    [
        ({"confidence": "high", "isHeterogeneous": "true"}, (None, None, None)),
        ({"confidence": ["High"], "isHeterogeneous": 1}, (None, None, None)),
        ({"confidence": 0.9, "isHeterogeneous": None}, (None, None, None)),
        (_review("Moderate", True, "  "), ("Moderate", True, None)),
        (_review("Moderate", True, 5), ("Moderate", True, None)),
        (_review("High", False, "mixed_lineage"), ("High", False, "mixed_lineage")),
        ({"confidence": "Low", "fullOutput": "n/a"}, ("Low", None, None)),
        ("not a review", (None, None, None)),
    ],
)
def test_reviewer_fields_are_type_checked(
    review: Any, expected: tuple[Any, Any, Any]
) -> None:
    latest = _new_result("1", _v2_cell_type())
    latest["review"] = review

    assert _sdk_values(_normalize(latest)) == dict(zip(SDK_KEYS, expected))


def test_missing_review_gives_none_values() -> None:
    latest = _new_result("1", _v2_cell_type())
    del latest["review"]

    assert _sdk_values(_normalize(latest)) == dict.fromkeys(SDK_KEYS)


@pytest.mark.parametrize(
    "cluster_data",
    [
        {},
        {"latest": {}},
        {"latest": {"developmentalSimilarity": []}},
        {"latest": {"annotation": "T cell", "confidence": 0.95}},
        {"latest": {"annotation": None}},
        {"latest": None},
        "not a cluster",
    ],
)
def test_clusters_without_an_annotation_are_skipped(cluster_data: Any) -> None:
    assert normalize_cluster_annotation("1", cluster_data) is None


@pytest.mark.parametrize(
    "latest",
    [
        _legacy_result("1"),
        _new_result("1", _gen2_cell_type()),
        _new_result("1", _v2_cell_type(), _review("Low", True, "doublet_artifact")),
        _new_result("1", _v2_cell_type(missingGenesReasoning="", conclusion="")),
    ],
)
def test_normalization_is_idempotent_and_does_not_modify_inputs(
    latest: dict[str, Any],
) -> None:
    cluster_data = {"latest": latest, "history": []}
    before = copy.deepcopy(cluster_data)

    entry = normalize_cluster_annotation("1", cluster_data)
    assert entry is not None
    entry_before = copy.deepcopy(entry)
    again = fill_cluster_annotation(entry, cluster_data)

    assert again == entry
    assert entry == entry_before
    assert cluster_data == before
    cell_type = (latest["annotation"]["fullOutput"] or {}).get("cellType", {})
    if "keySupportingGenes" in cell_type:
        # Filled lists are copies, so editing the result cannot change the raw data.
        assert entry["supportingMarkers"] is not cell_type["keySupportingGenes"]


def test_renormalize_results_fills_results_stored_by_older_versions() -> None:
    stored = _old_sdk_result(_results_response())
    before = copy.deepcopy(stored)

    result = renormalize_results(stored)
    annotations = _by_cluster(result)

    assert stored == before
    assert result is not stored
    assert _legacy_values(annotations["1"]) == LEGACY_EVIDENCE
    assert _legacy_values(annotations["2"]) == EXPECTED_FILL
    assert _legacy_values(annotations["3"]) == EXPECTED_FILL
    assert _sdk_values(annotations["3"]) == {
        "reviewerConfidence": "Low",
        "isHeterogeneous": True,
        "heterogeneityType": "distinct_populations",
    }
    assert result["raw_annotations"] == before["raw_annotations"]
    assert renormalize_results(result) == result


def test_renormalize_results_keeps_values_and_handles_unmatched_entries() -> None:
    stored = _old_sdk_result(_results_response())
    stored["annotations"][1]["justification"] = "Kept as stored"
    stored["annotations"][1]["reviewerConfidence"] = "High"
    stored["annotations"].append({"clusterId": "99", "annotation": "Orphan"})
    stored["annotations"].append("not an entry")

    result = renormalize_results(stored)

    assert result["annotations"][1]["justification"] == "Kept as stored"
    assert result["annotations"][1]["reviewerConfidence"] == "High"
    assert result["annotations"][1]["supportingMarkers"] == ["CD8A", "CD8B", "GZMK"]
    assert result["annotations"][3] == {
        "clusterId": "99",
        "annotation": "Orphan",
        **dict.fromkeys(SDK_KEYS),
    }
    assert result["annotations"][4] == "not an entry"


def test_renormalize_results_without_raw_annotations_adds_empty_review_fields() -> None:
    # Results saved by SDK 0.7.0 and older have no raw_annotations.
    stored = {
        "annotations": [
            {"clusterId": "1", "annotation": "T cell"},
            {"clusterId": "2", "annotation": "B cell", "reviewerConfidence": "High"},
            "not an entry",
        ]
    }
    before = copy.deepcopy(stored)

    result = renormalize_results(stored)

    assert result["annotations"] == [
        {"clusterId": "1", "annotation": "T cell", **dict.fromkeys(SDK_KEYS)},
        {
            "clusterId": "2",
            "annotation": "B cell",
            "reviewerConfidence": "High",
            "isHeterogeneous": None,
            "heterogeneityType": None,
        },
        "not an entry",
    ]
    assert result is not stored
    assert stored == before


def test_renormalize_results_pairs_duplicate_cluster_ids_in_order(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    response = {
        "studyContext": "",
        "annotations": {
            "a": {
                "latest": _new_result(
                    "X",
                    _v2_cell_type(keySupportingGenes=["CD8A"]),
                    _review("High", False),
                ),
                "history": [],
            },
            "no-result": {"latest": {}, "history": []},
            "b": {
                "latest": _new_result(
                    "X",
                    _v2_cell_type(keySupportingGenes=["MS4A1"]),
                    _review("Low", True, "mixed_lineage"),
                ),
                "history": [],
            },
        },
        "summary": {},
        "clusterCategories": [],
    }
    stored = _old_sdk_result(response)
    assert [entry["clusterId"] for entry in stored["annotations"]] == ["X", "X"]

    result = renormalize_results(stored)

    assert [entry["supportingMarkers"] for entry in result["annotations"]] == [
        ["CD8A"],
        ["MS4A1"],
    ]
    assert [entry["reviewerConfidence"] for entry in result["annotations"]] == [
        "High",
        "Low",
    ]
    # The same entries fetch_job_results builds from the response.
    assert result["annotations"] == _fetch(monkeypatch, response)["annotations"]


def test_renormalize_results_falls_back_to_the_raw_cluster_key() -> None:
    stored = {
        "annotations": [{"clusterId": "x", "annotation": "T cell"}],
        "raw_annotations": {
            "x": {"latest": {"annotation": "T cell", "review": {"confidence": "Low"}}}
        },
    }

    result = renormalize_results(stored)

    assert result["annotations"][0]["reviewerConfidence"] == "Low"


def test_load_local_results_renormalizes_h5ad_saved_by_older_versions(
    tmp_path: Path,
) -> None:
    old_result = _old_sdk_result(_results_response())
    adata = anndata.AnnData(obs=pd.DataFrame(index=["cell1", "cell2"]))
    adata.uns["cytetype_results"] = {
        "job_id": "old_job",
        "result": json.dumps(old_result),
    }
    path = tmp_path / "old.h5ad"
    adata.write_h5ad(path)
    adata = anndata.read_h5ad(path)
    stored_json = adata.uns["cytetype_results"]["result"]
    obs_columns = list(adata.obs.columns)

    result = load_local_results(adata, "cytetype")

    assert result is not None
    annotations = _by_cluster(result)
    assert _legacy_values(annotations["2"]) == EXPECTED_FILL
    assert _sdk_values(annotations["2"]) == {
        "reviewerConfidence": "Moderate",
        "isHeterogeneous": False,
        "heterogeneityType": None,
    }
    assert result["raw_annotations"] == old_result["raw_annotations"]
    assert adata.uns["cytetype_results"]["result"] == stored_json
    # Only the returned dict changes; no obs columns are backfilled.
    assert list(adata.obs.columns) == obs_columns


def test_load_local_results_does_not_mutate_dict_format_results() -> None:
    stored = _old_sdk_result(_results_response())
    before = copy.deepcopy(stored)
    adata = anndata.AnnData(obs=pd.DataFrame(index=["cell1"]))
    adata.uns["cytetype_results"] = {"job_id": "old_job", "result": stored}

    result = load_local_results(adata, "cytetype")

    assert result is not None
    assert result is not stored
    assert adata.uns["cytetype_results"]["result"] is stored
    assert stored == before
    assert result["annotations"][1] is not stored["annotations"][1]
    assert _legacy_values(result["annotations"][1]) == EXPECTED_FILL
    assert result["annotations"][1]["reviewerConfidence"] == "Moderate"


def test_store_annotations_adds_reviewer_confidence_column(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    log_messages: list[tuple[str, str]],
) -> None:
    result = _fetch(monkeypatch, _results_response())
    clusters = ["1", "1", "2", "2", "3", "4", "5", "3"]
    adata = anndata.AnnData(obs=pd.DataFrame(index=[f"cell{i}" for i in range(8)]))

    store_annotations(adata, result, "job-1", "cytetype", "leiden", clusters)
    messages = [message for _, message in log_messages]

    column = adata.obs["cytetype_reviewerConfidence_leiden"]
    assert isinstance(column.dtype, pd.CategoricalDtype)
    assert list(column.cat.categories) == ["High", "Moderate", "Low", "Unknown"]
    assert column.tolist() == [
        "High",
        "High",
        "Moderate",
        "Moderate",
        "Low",
        "Unknown",
        "Unknown",
        "Low",
    ]
    # The existing four columns keep their names and values.
    assert adata.obs["cytetype_annotation_leiden"].tolist() == [
        "Astrocyte",
        "Astrocyte",
        "CD8+ T Cell",
        "CD8+ T Cell",
        "CD8+ T Cell",
        "Unknown",
        "Unknown",
        "CD8+ T Cell",
    ]
    assert adata.obs["cytetype_cellOntologyTerm_leiden"].tolist()[:3] == [
        "astrocyte",
        "astrocyte",
        "CD8-positive, alpha-beta T cell",
    ]
    assert adata.obs["cytetype_cellOntologyTermID_leiden"].tolist()[-3:] == [
        "Unknown",
        "Unknown",
        "CL:0000625",
    ]
    assert adata.obs["cytetype_cellState_leiden"].tolist()[4:6] == ["Cytotoxic", ""]
    assert sorted(adata.obs.columns) == sorted(
        f"cytetype_{field}_leiden"
        for field in (
            "annotation",
            "cellOntologyTerm",
            "cellOntologyTermID",
            "cellState",
            "reviewerConfidence",
        )
    )
    summary = next(m for m in messages if "Full results" in m)
    assert "adata.obs['cytetype_cellOntologyTermID_leiden']" in summary
    assert "adata.obs['cytetype_reviewerConfidence_leiden']" in summary
    assert any("marked as 'Unknown'." in m for m in messages)

    path = tmp_path / "annotated.h5ad"
    adata.write_h5ad(path)
    reread = anndata.read_h5ad(path).obs["cytetype_reviewerConfidence_leiden"]
    assert list(reread.cat.categories) == ["High", "Moderate", "Low", "Unknown"]


def test_reviewer_confidence_column_keeps_only_occurring_categories() -> None:
    result = {
        "annotations": [
            {"clusterId": "1", "annotation": "T cell", "reviewerConfidence": "Low"},
            {"clusterId": "2", "annotation": "B cell", "reviewerConfidence": "bogus"},
            # Results stored without the key (older versions) read as Unknown.
            {"clusterId": "3", "annotation": "NK cell"},
        ]
    }
    adata = anndata.AnnData(obs=pd.DataFrame(index=["a", "b", "c"]))

    store_annotations(adata, result, "job-1", "cytetype", "leiden", ["1", "2", "3"])

    column = adata.obs["cytetype_reviewerConfidence_leiden"]
    assert list(column.cat.categories) == ["Low", "Unknown"]
    assert column.tolist() == ["Low", "Unknown", "Unknown"]


def _dotplot_response() -> dict[str, Any]:
    refined_v2 = _new_result(
        "0",
        _v2_cell_type(keySupportingGenes=["GENE1", "GENE2", "GENE3", "GENE4"]),
        label="Cytotoxic T Cell",
    )
    legacy = _legacy_result(
        "1", supportingMarkers=["GENE5", "GENE6"], annotation="Astrocyte"
    )
    del legacy["annotation"]["fullOutput"]
    return {
        "studyContext": "",
        "annotations": {
            # The server can replace annotation.annotation with a refined name.
            "0": {"latest": refined_v2, "history": []},
            "1": {"latest": legacy, "history": []},
            "2": {"latest": {"developmentalSimilarity": []}, "history": []},
        },
        "summary": {},
        "clusterCategories": [
            {"categoryName": "Lymphoid", "clusterIds": ["0"]},
            {"categoryName": "Glia", "clusterIds": ["1"]},
            {"categoryName": "Other", "clusterIds": ["2"]},
        ],
    }


def test_marker_dotplot_handles_legacy_results_and_unknown_clusters(
    mock_adata: anndata.AnnData,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    result = _fetch(monkeypatch, _dotplot_response())
    clusters = mock_adata.obs["leiden"].astype(str).tolist()
    store_annotations(mock_adata, result, "job-1", "cytetype", "leiden", clusters)

    dotplot = marker_dotplot(
        mock_adata,
        "leiden",
        gene_symbols="gene_symbols",
        show=False,
        return_fig=True,
    )
    dotplot.make_figure()
    group_labels = [text.get_text() for text in dotplot.ax_dict["gene_group_ax"].texts]
    plt.close("all")

    # "Unknown" has no category, so it goes last instead of aborting the plot.
    assert list(dotplot.categories_order) == [
        "Cytotoxic T Cell",
        "Astrocyte",
        "Unknown",
    ]
    assert list(dotplot.var_names) == ["GENE1", "GENE2", "GENE3", "GENE5", "GENE6"]
    # "Other" has no genes, so it gets no bracket or label.
    assert list(dotplot.var_groups.labels) == ["Lymphoid", "Glia"]
    assert list(dotplot.var_groups.positions) == [(0, 2), (3, 4)]
    assert group_labels == ["Lymphoid", "Glia"]


def test_marker_dotplot_requires_the_annotation_column(
    mock_adata: anndata.AnnData,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    result = _fetch(monkeypatch, _dotplot_response())
    mock_adata.uns["cytetype_results"] = {
        "job_id": "job-1",
        "result": json.dumps(result),
    }

    with pytest.raises(KeyError, match="cytetype_annotation_leiden"):
        marker_dotplot(mock_adata, "leiden")


def _dotplot_adata(
    monkeypatch: pytest.MonkeyPatch,
    var: pd.DataFrame,
    genes: dict[str, list[str]],
) -> anndata.AnnData:
    """AnnData annotated with clusters "0" (T cells) and "1" (B cells) and ``genes``."""
    labels = {"0": "T Cell", "1": "B Cell"}
    obs = pd.DataFrame(
        {"leiden": pd.Categorical(["0", "1"] * 10)},
        index=[f"cell{i}" for i in range(20)],
    )
    rng = np.random.default_rng(0)
    adata = anndata.AnnData(X=rng.random((20, len(var))), obs=obs, var=var)
    response = {
        "studyContext": "",
        "annotations": {
            cluster_id: {
                "latest": _new_result(
                    cluster_id,
                    _v2_cell_type(keySupportingGenes=cluster_genes),
                    label=labels[cluster_id],
                ),
                "history": [],
            }
            for cluster_id, cluster_genes in genes.items()
        },
        "summary": {},
        "clusterCategories": [
            {"categoryName": "T cells", "clusterIds": ["0"]},
            {"categoryName": "B cells", "clusterIds": ["1"]},
        ],
    }
    result = _fetch(monkeypatch, response)
    clusters = adata.obs["leiden"].astype(str).tolist()
    store_annotations(adata, result, "job-1", "cytetype", "leiden", clusters)
    return adata


def test_marker_dotplot_matches_mouse_genes_case_insensitively(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    var = pd.DataFrame(index=["Cd8a", "Gzmk", "Ms4a1", "Cd79a"])
    adata = _dotplot_adata(
        monkeypatch, var, {"0": ["CD8A", "GZMK"], "1": ["MS4A1", "CD79A"]}
    )

    dotplot = marker_dotplot(adata, "leiden", show=False, return_fig=True)

    assert list(dotplot.var_names) == ["Cd8a", "Gzmk", "Cd79a", "Ms4a1"]
    assert list(dotplot.var_groups.labels) == ["T cells", "B cells"]


def test_marker_dotplot_matches_the_gene_symbols_column(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    var = pd.DataFrame(
        {"symbol": ["Cd8a", "Gzmk", "Ms4a1", "Cd79a"]},
        index=[f"ENSMUSG{i:011d}" for i in range(4)],
    )
    adata = _dotplot_adata(
        monkeypatch, var, {"0": ["CD8A", "GZMK"], "1": ["MS4A1", "CD79A"]}
    )

    dotplot = marker_dotplot(
        adata, "leiden", gene_symbols="symbol", show=False, return_fig=True
    )

    assert list(dotplot.var_names) == ["Cd8a", "Gzmk", "Cd79a", "Ms4a1"]


def test_marker_dotplot_skips_unmatched_genes_with_one_warning(
    monkeypatch: pytest.MonkeyPatch,
    log_messages: list[tuple[str, str]],
) -> None:
    var = pd.DataFrame(index=["Cd8a", "Gzmk", "Ms4a1", "Cd79a"])
    adata = _dotplot_adata(
        monkeypatch,
        var,
        {"0": ["CD8A", "NOTAGENE"], "1": ["XIST1", "YGENE", "NOTAGENE"]},
    )

    dotplot = marker_dotplot(adata, "leiden", show=False, return_fig=True)

    assert list(dotplot.var_names) == ["Cd8a"]
    # B cells lost all its genes, so the group is dropped.
    assert list(dotplot.var_groups.labels) == ["T cells"]
    warnings = [message for level, message in log_messages if level == "WARNING"]
    assert warnings == [
        "3 marker genes are not shown because they were not found in "
        "adata.var_names: NOTAGENE, XIST1, YGENE."
    ]


def test_marker_dotplot_raises_when_no_marker_gene_is_found(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    var = pd.DataFrame(index=["Cd8a", "Gzmk"])
    adata = _dotplot_adata(monkeypatch, var, {"0": ["NOTAGENE"], "1": ["YGENE"]})

    with pytest.raises(ValueError, match="None of the marker genes were found"):
        marker_dotplot(adata, "leiden")
