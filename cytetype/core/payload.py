from typing import Any
import json
from pydantic import ValidationError

from ..api import LLMModelConfig, InputData, UploadedFiles
from ..config import logger


def build_annotation_payload(
    study_context: str,
    metadata: dict[str, Any] | None,
    cluster_map: dict[str, str],
    group_metadata: dict[str, dict[str, dict[str, int]]],
    marker_genes: dict[str, list[str]],
    visualization_data: dict[str, Any],
    expression_percentages: dict[str, dict[str, float]],
    n_parallel_clusters: int,
    llm_configs: list[dict[str, Any]] | None = None,
    uploaded_files: dict[str, str | None] | UploadedFiles | None = None,
) -> dict[str, Any]:
    """Build and validate API payload.

    Args:
        study_context: Biological context for the experimental setup
        metadata: Custom metadata tags for report header
        cluster_map: Mapping of cluster IDs to labels
        group_metadata: Aggregated metadata per cluster
        marker_genes: Marker genes per cluster
        visualization_data: Coordinates and cluster assignments
        expression_percentages: Gene expression percentages
        n_parallel_clusters: Number of parallel requests
        llm_configs: Optional LLM configurations
        uploaded_files: Optional upload IDs for obs_duckdb / vars_h5

    Returns:
        Validated payload dict for /annotate

    Raises:
        ValidationError: If payload validation fails
    """
    # Build and validate InputData directly
    try:
        validated_input = InputData(
            studyInfo=study_context,
            infoTags=metadata or {},
            clusterLabels={v: k for k, v in cluster_map.items()},
            clusterMetadata=group_metadata,
            markerGenes=marker_genes,
            visualizationData=visualization_data,
            expressionData=expression_percentages,
            nParallelClusters=n_parallel_clusters,
        )
        input_data = validated_input.model_dump()
    except ValidationError:
        logger.error(
            "Input data validation failed. This is likely due to incorrect data format "
            "from preprocessing. Please check that your AnnData object is properly formatted."
        )
        raise

    # Validate llm_configs
    validated_llm_configs: list[LLMModelConfig] = []
    if llm_configs:
        try:
            validated_llm_configs = [LLMModelConfig(**config) for config in llm_configs]
        except ValidationError:
            logger.error(
                "LLM configuration validation failed. Please check your llm_configs parameter. "
                "Each config must have 'provider' and 'name' fields, and either 'apiKey' or "
                "all AWS credentials ('awsAccessKeyId', 'awsSecretAccessKey', 'awsDefaultRegion')."
            )
            raise

    validated_uploaded_files: UploadedFiles | None = None
    if uploaded_files is not None:
        try:
            if isinstance(uploaded_files, UploadedFiles):
                validated_uploaded_files = uploaded_files
            else:
                validated_uploaded_files = UploadedFiles(**uploaded_files)
        except ValidationError:
            logger.error(
                "Uploaded file reference validation failed. Expected one or both "
                "of: obs_duckdb, vars_h5."
            )
            raise

    payload: dict[str, Any] = {
        "input_data": input_data,
        "llm_configs": (
            [config.model_dump() for config in validated_llm_configs]
            if validated_llm_configs
            else None
        ),
    }
    if validated_uploaded_files:
        payload["uploaded_files"] = validated_uploaded_files.model_dump(
            exclude_none=True
        )
    return payload


def save_query_to_file(input_data: dict[str, Any], filename: str) -> None:
    """Save query to JSON file.

    Args:
        input_data: The input data dictionary to save
        filename: Path to save the query file
    """
    with open(filename, "w") as f:
        json.dump(input_data, f)
