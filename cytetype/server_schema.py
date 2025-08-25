from typing import Any, TypeAlias, Literal
from pydantic import BaseModel, Field, model_validator


LLMProvider: TypeAlias = Literal[
    "google", "openai", "anthropic", "groq", "mistral", "openrouter", "bedrock"
]
AgentType: TypeAlias = Literal[
    "contextualizer", "annotator", "reviewer", "summarizer", "clinician", "chat"
]


class LLMModelConfig(BaseModel):
    provider: LLMProvider = Field(description="The provider of the model to use")
    name: str = Field(description="The name of the model to use")
    apiKey: str | None = Field(default=None, description="The API key for the model")
    baseUrl: str | None = Field(default=None, description="The base URL for the model")
    awsAccessKeyId: str | None = Field(
        default=None, description="The AWS access key ID for the model"
    )
    awsSecretAccessKey: str | None = Field(
        default=None, description="The AWS secret access key for the model"
    )
    awsDefaultRegion: str | None = Field(
        default=None, description="The AWS region for the model"
    )
    modelSettings: dict[str, Any] | None = Field(
        default=None, description="Extra body for the model"
    )
    targetAgents: list[AgentType] | None = Field(
        default_factory=list[AgentType],
        description="List of agents that can use this model",
    )
    skipValidation: bool = Field(
        default=False,
        description="Whether to skip validation for the model. Only turn this off if you are sure that the model has worked before.",
    )
    allowFallback: bool = Field(
        default=False,
        description="Whether to allow the model to fallback to a default model",
    )

    # add a model validator to check that all the aws credentials are provided if any of them are provided. also check that either all aws credentials are provided or none of them are provided. or apikey is provided.
    @model_validator(mode="after")
    def check_aws_credentials(self) -> "LLMModelConfig":
        if self.apiKey is not None:
            return self
        if self.awsAccessKeyId is not None:
            if self.awsSecretAccessKey is None or self.awsDefaultRegion is None:
                raise ValueError(
                    "All AWS credentials must be provided if any of them are provided"
                )
            return self
        raise ValueError("Either apiKey or all AWS credentials must be provided")


class InputData(BaseModel):
    studyInfo: str = Field(
        default="",
        description="The biological context for the experimental setup. Can include information about organisms, tissues, diseases, developmental stages, single-cell methods, experimental conditions",
    )
    infoTags: dict[str, str] = Field(
        default_factory=dict,
        description="Optional metadata key-value pairs to display in the report header. Values that look like URLs will be made clickable.",
    )
    clusterLabels: dict[str, str] = Field(
        default_factory=dict,
        description="Optional user-provided labels for cluster IDs to improve navigation. Keys are cluster IDs, values are display labels.",
    )
    clusterMetadata: dict[str, dict[str, dict[str, int]]] = Field(
        default_factory=dict,
        description="This is a dictionary of dictionaries of dictionaries. The first key is the cluster id, the second key is the one of metadata columns, and the third key is the percentage of cells in the cluster that have the metadata value.",
    )
    markerGenes: dict[str, list[str]] = Field(
        default_factory=dict,
        description="Dictionary mapping cluster IDs to their marker genes",
    )
    visualizationData: dict[str, Any] | None = Field(
        default=None,
        description="Optional visualization data containing coordinates and cluster assignments for scatter plots",
    )
    expressionData: dict[str, dict[str, float]] = Field(
        default_factory=dict,
        description="Dictionary mapping gene names to their expression percentages across clusters",
    )

    @classmethod
    def get_example(cls) -> "InputData":
        """Return an example InputData instance for documentation and testing"""
        return cls(
            studyInfo="Adult human brain tissue samples from healthy controls and Alzheimer's disease patients, analyzed using 10X Genomics single-cell RNA-seq. Samples include cortical and hippocampal regions.",
            infoTags={
                "Study": "Alzheimer's Disease Brain Atlas",
                "DOI": "https://doi.org/10.1038/s41586-023-06063-y",
                "GEO Accession": "GSE157827",
                "PubMed": "https://pubmed.ncbi.nlm.nih.gov/37258686/",
                "Dataset": "10X Genomics scRNA-seq",
                "Tissue": "Human brain cortex and hippocampus",
            },
            clusterLabels={
                "Cluster1": "Astrocytes",
                "Cluster2": "Neurons",
            },
            clusterMetadata={
                "Cluster1": {
                    "condition": {"healthy": 60, "alzheimer": 40},
                    "region": {"cortex": 70, "hippocampus": 30},
                    "donor": {"donor1": 25, "donor2": 35, "donor3": 40},
                },
                "Cluster2": {
                    "condition": {"healthy": 80, "alzheimer": 20},
                    "region": {"cortex": 90, "hippocampus": 10},
                    "donor": {"donor1": 30, "donor2": 30, "donor3": 40},
                },
            },
            markerGenes={
                "Cluster1": ["GFAP", "S100B", "AQP4"],
                "Cluster2": ["RBFOX3", "MAP2", "SYP"],
            },
            visualizationData={
                "coordinates": [
                    [-2.1, 1.5],
                    [3.2, -0.8],
                    [1.7, 2.3],
                    [-1.1, -1.9],
                    [2.8, 1.1],
                    [-0.5, 3.1],
                    [4.2, -2.1],
                    [0.3, 0.7],
                ],
                "clusters": [
                    "Cluster1",
                    "Cluster1",
                    "Cluster2",
                    "Cluster2",
                    "Cluster1",
                    "Cluster2",
                    "Cluster1",
                    "Cluster2",
                ],
            },
            expressionData={
                "GFAP": {
                    "Cluster1": 85.2,
                    "Cluster2": 3.1,
                    "Cluster3": 4.5,
                },
                "S100B": {
                    "Cluster1": 92.7,
                    "Cluster2": 2.8,
                    "Cluster3": 5.2,
                },
            },
        )
