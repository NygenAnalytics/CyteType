# Configuration

## Initialization Parameters
```python
from cytetype import CyteType

annotator = CyteType(
    adata,
    group_key="leiden",                 # Required: cluster column name
    rank_key="rank_genes_groups",       # DE results key
    gene_symbols_column="gene_symbols", # Gene symbols column
    n_top_genes=50,                      # Top marker genes per cluster
    aggregate_metadata=True,             # Aggregate metadata
    min_percentage=10,                   # Min percentage for cluster context
    pcent_batch_size=2000,               # Batch size for calculations
    coordinates_key="X_umap",           # Coordinates key for visualization
    max_cells_per_group=1000,            # Max cells per group for visualization
)
```

## Run Options
```python
adata = annotator.run(
    study_context="Human PBMC from COVID-19 patients",
    metadata={"experiment_name": "PBMC_COVID", "run_label": "v1"},
    save_query=True,
    query_filename="query.json",
    show_progress=True,
)
```

## LLM Configuration
You can provide your own LLM providers/models:
```python
adata = annotator.run(
    study_context="Human PBMC",
    llm_configs=[{
        "provider": "openai",   # one of: anthropic, bedrock, google, groq, mistral, openai, openrouter
        "name": "gpt-4o-mini",
        "apiKey": "your-api-key",
        "baseUrl": "https://api.openai.com/v1",   # optional
        "modelSettings": {                        # optional
            "temperature": 0.0,
            "max_tokens": 4096
        }
    }],
)
```
- For local models via Ollama, see [Ollama Integration](./ollama.md)

## Authentication

Before the first annotation, sign in through the CLI:

```bash
cytetype setup
```

CyteType opens browser-based authorization and saves a personal API key locally. The Python client automatically loads the saved key for the selected API server when `run()` starts.

For a custom server, use the same origin during setup and in Python:

```bash
cytetype setup --api-url https://cytetype.example.org
```

```python
annotator = CyteType(
    adata,
    group_key="leiden",
    api_url="https://cytetype.example.org",
)
adata = annotator.run(study_context="Human PBMC")
```

For CI or another managed environment, a token can be supplied directly instead of using the local credentials file:

```python
import os

adata = annotator.run(
    study_context="Human PBMC",
    auth_token=os.environ["CYTETYPE_API_TOKEN"],
)
```

The environment variable in this example is user-managed. CyteType does not read `CYTETYPE_API_TOKEN` automatically. Do not hard-code API keys in source files or notebooks.

See [CLI and Authentication](./cli.md) for all commands, credential locations, custom server rules, and authentication precedence.

## Artifacts

`run()` automatically builds and uploads two artifact files before submitting an annotation job:

- **`vars.h5`**: a compressed HDF5 file containing the normalized expression matrix (`adata.X`) and variable metadata (`adata.var`). Used by the server for on-demand gene expression lookups during annotation and in the interactive report.
- **`obs.duckdb`**: a DuckDB database containing the observation metadata (`adata.obs`). Used by the server to power metadata queries and filtering in the interactive report.

Both files are created locally and then uploaded to the CyteType API. The uploaded references are attached to the `/annotate` payload so the server can link them to the job.

### Artifact Parameters

```python
adata = annotator.run(
    ...
    vars_h5_path="vars.h5",        # Local output path for vars artifact
    obs_duckdb_path="obs.duckdb",   # Local output path for obs artifact
    upload_timeout_seconds=3600,    # Socket read timeout per upload (seconds)
    cleanup_artifacts=False,        # Delete local artifact files after run
    require_artifacts=True,         # Raise on artifact failure (set False to skip)
)
```

| Parameter | Default | Description |
|-----------|---------|-------------|
| `vars_h5_path` | `"vars.h5"` | Local path where the vars HDF5 file is written |
| `obs_duckdb_path` | `"obs.duckdb"` | Local path where the obs DuckDB file is written |
| `upload_timeout_seconds` | `3600` | Socket read timeout for each artifact upload |
| `cleanup_artifacts` | `False` | Delete local artifact files after run completes or fails |
| `require_artifacts` | `True` | Raise on artifact build/upload failure. Set to `False` to skip artifacts and continue with annotation only |

### Error Handling

By default (`require_artifacts=True`), any failure during artifact building or uploading stops the run and surfaces the full error. The error message includes a link to report the issue on GitHub.

If you want the annotation to proceed even when artifacts fail (e.g. due to disk space or network issues), set `require_artifacts=False`. Annotation will still work, but the interactive report will not have expression lookups or metadata filtering.

### Memory Recommendation for Large Datasets

For large datasets, open your AnnData object in backed mode to reduce memory usage while building `vars.h5`:

```python
import scanpy as sc

adata = sc.read_h5ad("input.h5ad", backed="r")
```

## Advanced

```python
adata = annotator.run(
    ...
    poll_interval_seconds=30,   # How often to poll (default)
    timeout_seconds=7200,       # Max wait time (default: 2 hours)
    api_url="https://custom-api.example.com",  # Custom API endpoint if needed
)
```
