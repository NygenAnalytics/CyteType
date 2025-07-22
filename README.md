<h1 align="left">CyteType</h1>

<p align="left">
  <!-- GitHub Actions CI Badge -->
  <a href="https://github.com/NygenAnalytics/cytetype/actions/workflows/publish.yml">
    <img src="https://github.com/NygenAnalytics/cytetype/actions/workflows/publish.yml/badge.svg" alt="CI Status">
  </a>
  <a href="https://github.com/NygenAnalytics/cytetype/blob/main/LICENSE">
    <img src="https://img.shields.io/badge/License-CC%20BY--NC--SA%204.0-lightgrey.svg" alt="License: CC BY-NC-SA 4.0">
  </a>
  <a href="https://pypi.org/project/cytetype/">
    <img src="https://img.shields.io/pypi/v/cytetype.svg" alt="PyPI version">
  </a>
  <img src="https://img.shields.io/badge/python-≥3.11-blue.svg" alt="Python Version">
  <a href="https://colab.research.google.com/drive/1aRLsI3mx8JR8u5BKHs48YUbLsqRsh2N7?usp=sharing">
    <img src="https://colab.research.google.com/assets/colab-badge.svg" alt="Open In Colab">
  </a>
</p>

---

> ⚠️ CyteType is under active development and breaking changes may be introduced. Please work with the latest version to ensure compatibility and access to new features.

**CyteType** is a Python package for deep characterization of cell clusters from single-cell RNA-seq data. This package interfaces with Anndata objects to call CyteType API.

## Table of Contents

- [Example Report](#example-report)
- [Quick Start](#quick-start)
- [Installation](#installation)
- [Usage](#usage)
  - [Required Preprocessing](#required-preprocessing)
  - [Annotation](#annotation)
- [Configuration Options](#configuration-options)
  - [Initialization Parameters](#initialization-parameters)
  - [Submitting Annotation job](#submitting-annotation-job)
  - [Custom LLM Configuration](#custom-llm-configuration)
  - [Custom LLM Configuration (Ollama)](#custom-llm-configuration-ollama)
  - [Advanced parameters](#advanced-parameters)
- [Annotation Process](#annotation-process)
  - [Core Functionality](#core-functionality)
  - [Advanced Context Generation](#advanced-context-generation)
  - [Result Format](#result-format)
  - [Advanced Result Components](#advanced-result-components)
- [Development](#development)
  - [Setup](#setup)
  - [Exception Handling](#exception-handling)
  - [Testing](#testing)
- [License](#license)

## Example Report

View a sample annotation report: <a href="https://nygen-labs-prod--cytetype-api.modal.run/report/77069508-d9f1-4a79-bdab-5870fc3ccdf3?v=250722" target="blank">CyteType Report</a>



## Quick Start

```python
import anndata
import scanpy as sc
import cytetype

# Load and preprocess your data
adata = anndata.read_h5ad("path/to/your/data.h5ad")
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.pca(adata)
sc.pp.neighbors(adata)
sc.tl.leiden(adata, key_added = "clusters") 
sc.tl.rank_genes_groups(adata, groupby='clusters', method='t-test')

# Initialize CyteType (performs data preparation)
annotator = cytetype.CyteType(adata, group_key='clusters')

# Run annotation
adata = annotator.run(
    study_context="Human brain tissue from Alzheimer's disease patients"
)

# View results
print(adata.obs.cytetype_annotation_clusters)
print(adata.obs.cytetype_cellOntologyTerm_clusters)
```

## Installation

```bash
pip install cytetype
```

## Usage

### Required Preprocessing

Your `AnnData` object must have:

- Log-normalized expression data in `adata.X`
- Cluster labels in `adata.obs` 
- Differential expression results from `sc.tl.rank_genes_groups`

```python
import scanpy as sc

# Standard preprocessing
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# Clustering
sc.pp.pca(adata)
sc.pp.neighbors(adata)
sc.tl.leiden(adata, key_added='clusters')

# Differential expression (required)
sc.tl.rank_genes_groups(adata, groupby='clusters', method='t-test')
```

### Annotation

```python
from cytetype import CyteType

# Initialize (data preparation happens here)
annotator = CyteType(adata, group_key='clusters')

# Run annotation
adata = annotator.run(
    study_context="Adult human brain tissue samples from healthy controls and Alzheimer's disease patients, analyzed using 10X Genomics single-cell RNA-seq. Samples include cortical and hippocampal regions."
)

# Or with custom metadata for tracking
adata = annotator.run(
    study_context="Adult human brain tissue samples from healthy controls and Alzheimer's disease patients, analyzed using 10X Genomics single-cell RNA-seq. Samples include cortical and hippocampal regions.",
    metadata={
        'experiment_name': 'Brain_AD_Study',
        'run_label': 'initial_analysis'
    }
)

# Results are stored in:
# - adata.obs.cytetype_annotation_clusters (cell type annotations)
# - adata.obs.cytetype_cellOntologyTerm_clusters (cell ontology terms)
# - adata.uns['cytetype_results'] (full API response)
```

# Configuration Options

## Initialization Parameters

```python
annotator = CyteType(
    adata,
    group_key='leiden',                    # Required: cluster column name
    rank_key='rank_genes_groups',          # DE results key
    gene_symbols_column='gene_symbols',    # Gene symbols column
    n_top_genes=50,                        # Top marker genes per cluster
    aggregate_metadata=True,               # Aggregate metadata
    min_percentage=10,                     # Min percentage for cluster context
    pcent_batch_size=2000,                 # Batch size for calculations
    coordinates_key='X_umap',              # Coordinates key for visualization
    max_cells_per_group=1000,              # Max cells per group for visualization 
)
```

## Submitting Annotation job

The `run` method accepts several configuration parameters to control the annotation process:

```python
annotator.run(
    study_context="Adult human brain tissue samples from healthy controls and Alzheimer's disease patients, analyzed using 10X Genomics single-cell RNA-seq. Samples include cortical and hippocampal regions.",
    metadata={
        'experiment_name': 'Brain_AD_Study',
        'run_label': 'initial_analysis'
    },
    save_query=True,
    query_filename="query.json",
    show_progress=True,
)
```

### Custom LLM Configuration

The CyteType API provides access to some chosen LLM providers by default.
Users can choose to provide their own LLM models and model providers.
Many models can be provided simultaneously, and then they will be used iteratively for each of the clusters.

```python
adata = annotator.run(
    study_context="Human PBMC from COVID-19 patients",
    llm_configs=[{
        'provider': 'openai',
        'name': 'gpt-4o-mini',
        'apiKey': 'your-api-key',
        'baseUrl': 'https://api.openai.com/v1',  # Optional
        'modelSettings': {                       # Optional
            'temperature': 0.0,
            'max_tokens': 4096
        }  
    }],
)
```

#### Rate Limits

If you do not provide your own model providers, then the CyteType API implements rate limiting for fair usage:

- Annotation submissions: 5 RPD
- Reannotation: 10 RPD
- Report retrieval: 20 RPM

If you exceed rate limits, the system will return appropriate error messages with retry timing information

Supported providers: `openai`, `anthropic`, `google`, `xai`, `groq`, `mistral`, `openrouter`, `bedrock`

### Custom LLM Configuration (Ollama)

The CyteType API supports Ollama models as well. You will need to expose your Ollama server to the internet using a tunneling service. Refer to the [OLLAMA.md](./OLLAMA.md) file for instructions on how to do this.

### Advanced parameters

```python
adata = annotator.run(
    ...
    # API polling and timeout settings
    poll_interval_seconds=30,           # How often to check for results (default)
    timeout_seconds=7200,               # Max wait time (default: 2 hours)
    
    # API configuration
    api_url="https://custom-api.com",   # Custom API endpoint
    auth_token="your-auth-token",       # Authentication token
)
```

### Authentication and Authorization

You can provide your own token to the `run` method using the `auth_token` parameter. This will be included in the Authorization header as "Bearer {auth_token}". All API requests will be authenticated with this token.

## Annotation Process

CyteType performs comprehensive cell type annotation through an automated pipeline:

### Core Functionality

- **Automated Annotation**: Identifies likely cell types for each cluster based on marker genes
- **Ontology Mapping**: Maps identified cell types to Cell Ontology terms (e.g., `CL_0000127`)  
- **Review & Justification**: Analyzes supporting/conflicting markers and assesses confidence
- **Literature Search**: Searches for relevant literature to support the annotation

### Advanced Context Generation

CyteType generates detailed contextual information to inform annotations:

**Dataset-Level Context**: Comprehensive analysis of experimental metadata:
```
"This dataset originates from multiple human tissues including adrenal gland, 
brain, liver, lung, lymph node, and pleural effusion, with samples derived 
from both healthy individuals and patients diagnosed with lung adenocarcinoma 
or small cell lung carcinoma. Experimental data was generated via 10X Genomics 
Chromium 3' single-cell sequencing, which may introduce platform-specific 
technical artifacts."
```

**Cluster-Specific Context**: Detailed metadata analysis for each cluster:
```
"Cluster 1 comprises 99% lung-derived cells, with 65% originating from lung 
adenocarcinoma samples and 33% from normal tissue. The cells are distributed 
across two primary donors with demographic characteristics including 67% 
female donors and 97% self-reported European ethnicity. Treatment conditions 
include Platinum Doublet (55%) and Naive (44%)."
```

This contextual information enables more accurate annotations by considering:
- **Tissue Origins**: Multi-tissue datasets with precise anatomical mapping
- **Disease States**: Healthy vs. pathological conditions with treatment history
- **Technical Factors**: Sequencing platforms, batch effects, and processing methods
- **Demographics**: Age, sex, and ethnicity distributions
- **Treatment Context**: Therapeutic interventions and their potential cellular effects

### Result Format

Results include comprehensive annotations for each cluster with expert-level analysis:

```python
# Access results after annotation using the helper method
results = annotator.get_results()

# Or access directly from the stored JSON string
import json
results = json.loads(adata.uns['cytetype_results']['result'])

# Each annotation includes comprehensive information:
for annotation in results['annotations']:
    print(f"Cluster: {annotation['clusterId']}")
    print(f"Cell Type: {annotation['annotation']}")
    print(f"Granular Annotation: {annotation['granularAnnotation']}")
    print(f"Cell State: {annotation['cellState']}")
    print(f"Confidence: {annotation['confidence']}")
    print(f"Ontology Term: {annotation['ontologyTerm']}")
    print(f"Is Heterogeneous: {annotation['isHeterogeneous']}")
    
    # Supporting evidence and conflicts
    print(f"Supporting Markers: {annotation['supportingMarkers']}")
    print(f"Conflicting Markers: {annotation['conflictingMarkers']}")
    print(f"Missing Expression: {annotation['missingExpression']}")
    print(f"Unexpected Expression: {annotation['unexpectedExpression']}")
    
    # Expert review and justification
    print(f"Justification: {annotation['justification']}")
    print(f"Review Comments: {annotation['reviewComments']}")
    print(f"Feedback: {annotation['feedback']}")
    
    # Similarity and literature support
    print(f"Similar Clusters: {annotation['similarity']}")
    print(f"Corroborating Papers: {len(annotation['corroboratingPapers']['papers'])} papers")
    
    # Model usage and performance metrics
    print(f"Models Used: {annotation['llmModels']}")
    print(f"Total Processing Time: {annotation['usageInfo']['total_runtime_seconds']:.1f}s")
    print(f"Total Tokens: {annotation['usageInfo']['total_tokens']}")
```

#### Advanced Result Components

**Expert Review System**: Each annotation undergoes multi-stage review with detailed feedback:
- **Review Comments**: Expert-level biological interpretation and mechanistic insights
- **Confidence Assessment**: Moderate/High confidence based on marker evidence
- **Feedback Loop**: Iterative refinement based on biological plausibility
- **Mechanistic Analysis**: Discussion of signaling pathways, developmental biology, and disease pathogenesis

**Literature Integration**: Automatic literature search provides supporting evidence:
- **Corroborating Papers**: Relevant publications with PMIDs and summaries
- **Biological Context**: Integration of current research to validate annotations

Example corroborating papers:
```python
papers = annotation['corroboratingPapers']['papers']
for paper in papers:
    print(f"Title: {paper['title']}")
    print(f"PMID: {paper['pmid']}")
    print(f"Journal: {paper['journal']} ({paper['year']})")
    print(f"Summary: {paper['summary']}")
```

Sample output:
```
Title: YAP regulates alveolar epithelial cell differentiation and AGER via NFIB/KLF5/NKX2-1
PMID: 34466790
Journal: iScience (2021)
Summary: Documents atypical HOPX+AGER+SFTPC+ 'dual-positive' alveolar cells that 
persist in mature lungs, directly validating the mixed AT1/AT2 phenotype observed 
in malignant clusters.
```

**Marker Analysis**: Comprehensive evaluation of gene expression patterns:
- **Supporting Markers**: Genes that strongly support the annotation
- **Conflicting Markers**: Genes that challenge the annotation with explanations
- **Missing/Unexpected Expression**: Detailed analysis of expression anomalies with biological explanations

Example unexpected expression analysis:
```
"Expression of AT2 markers (SFTPC, SFTPB) and club cell marker (SCGB1A1) 
in a cluster with strong AT1 markers" 
→ Explained by: "dedifferentiation process in cancer where transformed 
epithelial cells exhibit aberrant co-expression of markers from multiple 
lineages due to pathological plasticity"
```

**Performance Metrics**: Detailed usage statistics for transparency:
- **Model Information**: Which LLM models were used for each analysis step
- **Runtime Statistics**: Processing time and token usage per cluster
- **Annotation Attempts**: Number of refinement iterations

#### Example Annotations

CyteType provides sophisticated, multi-layered annotations:

**Basic Cell Type**: `"B-cell"`, `"Lung Adenocarcinoma Cell"`

**Granular Annotations**: Detailed phenotypic descriptions:
- `"AGER-positive, HOPX-positive, KRT19-positive lung adenocarcinoma cell with mixed AT1/AT2 phenotype"`
- `"EMT-transitioned, pleural metastasis-competent adenocarcinoma cell with platinum-induced stress phenotype"`
- `"CD74-high activated tumor-infiltrating B-cell in lung adenocarcinoma microenvironment"`

**Cell States**: Functional and pathological states:
- `"Transformed"`, `"Malignant"`, `"Activated"`, `"EMT-transitioned and stressed"`

**Expert Review Comments**: Detailed mechanistic insights:
```
"The mixed AT1/AT2 phenotype observed in this malignant cluster exemplifies 
the pathological dedifferentiation characteristic of lung adenocarcinoma, 
but the degree of lineage promiscuity suggests exceptional cellular plasticity 
beyond typical adenocarcinoma patterns. This may indicate activation of 
primitive developmental pathways like Wnt/β-catenin signaling..."
```

## Development

### Setup

```bash
git clone https://github.com/NygenAnalytics/cytetype.git
cd cytetype
uv sync --all-extras
uv run pip install -e .
```

### Exception Handling

The package defines several custom exceptions for different error scenarios:

- **`CyteTypeError`**: Base exception class for all CyteType-related errors
- **`CyteTypeAPIError`**: Raised for errors during API communication (network issues, invalid responses)
- **`CyteTypeTimeoutError`**: Raised when API requests timeout
- **`CyteTypeJobError`**: Raised when the API reports an error for a specific job

### Testing

```bash
uv run pytest              # Run tests
uv run ruff check .        # Linting
uv run ruff format .       # Formatting
uv run mypy .              # Type checking
```

## License

Licensed under CC BY-NC-SA 4.0 - see [LICENSE](LICENSE) for details.
