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
```python
adata = annotator.run(
    study_context="...",
    auth_token="your-auth-token",  # included as Authorization: Bearer <token>
)
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