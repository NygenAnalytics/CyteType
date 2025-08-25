# Troubleshooting

- Rate limit responses include retry information; wait or provide your own LLM.
- Verify preprocessing: clustering and `rank_genes_groups` must be present. 
- Make sure you have valid gene symbols in the AnnData object and are passing the correct gene symbols column name to parameter `gene_symbols_column`.
- If you are using a custom LLM, make sure you have the correct API key and base URL.