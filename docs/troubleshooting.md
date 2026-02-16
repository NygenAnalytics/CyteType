# Troubleshooting

- Rate limit responses include retry information; wait or provide your own LLM.
- Verify preprocessing: clustering and `rank_genes_groups` must be present. 
- Make sure you have valid gene symbols in the AnnData object and are passing the correct gene symbols column name to parameter `gene_symbols_column`.
- If you are using a custom LLM, make sure you have the correct API key and base URL.
- For large datasets, load AnnData in backed mode (`sc.read_h5ad(..., backed="r")`) to reduce memory use during artifact generation.
- `run()` creates `vars.h5` and `obs.duckdb` before annotation. Use `cleanup_artifacts=True` if you do not want to keep these local files.