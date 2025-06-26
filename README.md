<h1 align="left">CyteType</h1>

<p align="left">
  <!-- GitHub Actions CI Badge -->
  <a href="https://github.com/NygenAnalytics/cytetype/actions/workflows/publish.yml">
    <img src="https://github.com/NygenAnalytics/cytetype/actions/workflows/publish.yml/badge.svg" alt="CI Status">
  </a>
  <a href="https://github.com/NygenAnalytics/cytetype/blob/main/LICENSE">
    <img src="https://img.shields.io/badge/License-CC%20BY--NC--SA%204.0-lightgrey.svg" alt="License: CC BY-NC-SA 4.0">
  </a>
  <img src="https://img.shields.io/badge/python-≥3.11-blue.svg" alt="Python Version">
  <a href="https://colab.research.google.com/drive/1aRLsI3mx8JR8u5BKHs48YUbLsqRsh2N7?usp=sharing">
    <img src="https://colab.research.google.com/assets/colab-badge.svg" alt="Open In Colab">
  </a>
</p>

---

> **⚠️ Important Notice - URL Update**: The CyteType API URL has been updated. If you have bookmarked or saved any old report links, please update them to use the new domain. The new API endpoint is now active at `https://cytetype.nygen.io`.

**CyteType** is a Python package for deep chracterization of cell clusters from single-cell RNA-seq data. This package interfaces with Anndata objects to call CyteType API.

## Example Report

View a sample annotation report: [CyteType Report](https://cytetype.nygen.io/report/6263d2ba-0865-4edb-bec1-06b35be4e80b)

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

The `study_context` should include comprehensive biological information about your experimental setup:

- **Organisms**: Species being studied (e.g., "human", "mouse")
- **Tissues**: Tissue types and anatomical regions
- **Diseases**: Disease conditions or states
- **Developmental stages**: Age, developmental timepoints
- **Single-cell methods**: Sequencing platform (e.g., "10X Genomics", "Smart-seq2")
- **Experimental conditions**: Treatments, time courses, perturbations

**Example**: `"Adult human brain tissue samples from healthy controls and Alzheimer's disease patients, analyzed using 10X Genomics single-cell RNA-seq. Samples include cortical and hippocampal regions."`

# Configuration Options

## Initialization Parameters

```python
annotator = CyteType(
    adata,
    group_key='leiden',                    # Required: cluster column name
    rank_key='rank_genes_groups',          # DE results key (default)
    gene_symbols_column='gene_symbols',    # Gene symbols column (default)
    n_top_genes=50,                        # Top marker genes per cluster
    aggregate_metadata=True,               # Aggregate metadata (default)
    min_percentage=10,                     # Min percentage for cluster context
    pcent_batch_size=2000,                 # Batch size for calculations
)
```

## Submitting Annotation job

The `run` method accepts several configuration parameters to control the annotation process:

### Custom LLM Configuration

The CyteType API provides access to some chosen LLM providers by default.
Users can choose to provide their own LLM models and model providers.
Many models can be provided simultaneously and then they will be used iteratively for each of the clusters.

```python
adata = annotator.run(
    study_context="Human PBMC from COVID-19 patients",
    model_config=[{
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

If you do not provide your own model providers then the CyteType API implements rate limiting for fair usage:

- Annotation submissions: 5 requests per hour per IP
- Result retrieval: 20 requests per minute per IP

If you exceed rate limits, the system will return appropriate error messages with retry timing information

Supported providers: `openai`, `anthropic`, `google`, `xai`, `groq`, `mistral`, `openrouter`

### Custom LLM Configuration (Ollama)

The CyteType API supports Ollama models as well. You will need to expose your Ollama server to the internet using a tunneling service. Refer to the [OLLAMA.md](./OLLAMA.md) file for instructions on how to do this.

### Advanced parameters

```python
adata = annotator.run(
    ...
    run_config={
        'concurrentClusters': 3,        # Default: 3, Range: 2-10
        'maxAnnotationRevisions': 2,    # Default: 2, Range: 1-5
    },
    
    # Custom metadata for tracking
    metadata={
        'experiment_name': 'PBMC_COVID_Study',
        'run_label': 'baseline_analysis',
        'researcher': 'Dr. Smith',
        'batch': 'batch_001'
    },
    
    # API polling and timeout settings
    poll_interval_seconds=10,           # How often to check for results
    timeout_seconds=1200,               # Max wait time (20 minutes)
    
    # API configuration
    api_url="https://custom-api.com",   # Custom API endpoint
    auth_token="your-auth-token",       # Authentication token
    save_query=True                     # Save query to query.json
)
```

#### Run configuration

- **`concurrentClusters`** (int, default=5, range=2-30): Maximum number of clusters to process simultaneously. Higher values may speed up processing but can cause rate limit errors from LLM API providers.
- **`maxAnnotationRevisions`** (int, default=2, range=1-5): Maximum number of refinement iterations based on reviewer feedback. More revisions may improve annotation quality but increase processing time.

#### Additional Run Parameters

- **`metadata`** (dict, optional): Custom metadata to send with the API request for tracking purposes. Can include experiment names, run labels, researcher information, or any other user-defined data. This metadata is sent to the API but not stored locally with results.
- **`poll_interval_seconds`** (int, default=30): How frequently to check the API for job completion.
- **`timeout_seconds`** (int, default=3600): Maximum time to wait for results before timing out.
- **`api_url`** (str): Custom API endpoint URL for self-hosted deployments.
- **`auth_token`** (str, optional): Bearer token for API authentication.
- **`save_query`** (bool, default=True): Whether to save the API query to `query.json` for debugging.

## Annotation Process

CyteType performs comprehensive cell type annotation through an automated pipeline:

### Core Functionality

- **Automated Annotation**: Identifies likely cell types for each cluster based on marker genes
- **Ontology Mapping**: Maps identified cell types to Cell Ontology terms (e.g., `CL_0000127`)  
- **Review & Justification**: Analyzes supporting/conflicting markers and assesses confidence
- **Alternative Suggestions**: Provides potential alternative annotations when applicable

### Result Format

Results include detailed annotations for each cluster:

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
    print(f"Is Approved: {annotation['is_approved']}")
    print(f"Is Heterogeneous: {annotation['isHeterogeneous']}")
    print(f"Supporting Markers: {annotation['supportingMarkers']}")
    print(f"Conflicting Markers: {annotation['conflictingMarkers']}")
    print(f"Missing Expression: {annotation['missingExpression']}")
    print(f"Unexpected Expression: {annotation['unexpectedExpression']}")
    print(f"Alternative Annotations: {annotation['alternativeAnnotations']}")
    print(f"Justification: {annotation['justification']}")
    print(f"Review Comments: {annotation['reviewComments']}")
    print(f"Feedback: {annotation['feedback']}")
    print(f"Similarity: {annotation['similarity']}")
```

## Development

### Setup

```bash
git clone https://github.com/NygenAnalytics/cytetype.git
cd cytetype
uv sync --all-extras
uv run pip install -e .
```

### Testing

```bash
uv run pytest              # Run tests
uv run ruff check .        # Linting
uv run ruff format .       # Formatting
uv run mypy .              # Type checking
```

## License

Licensed under CC BY-NC-SA 4.0 - see [LICENSE](LICENSE) for details.
