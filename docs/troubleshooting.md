# Troubleshooting

## CLI and Authentication

- If Python reports `CyteType sign-in is required`, run `cytetype setup`, complete browser authorization, and retry.
- If you already have an API key, run `cytetype login`. The key is entered through a hidden prompt and saved only after the server validates it.
- For a custom server, use the same origin in the CLI and Python. For example, run `cytetype setup --api-url https://cytetype.example.org` and pass `api_url="https://cytetype.example.org"` to `CyteType`.
- If the browser does not open during setup, copy the authorization URL printed in the terminal. Keep the command running while you authorize because it is waiting for a callback on `127.0.0.1`.
- If setup times out, retry and complete authorization within five minutes. Check whether a local firewall or browser policy is blocking the callback to `127.0.0.1`.
- If setup reports that the saved key is invalid or inactive, run `cytetype setup --force` to authenticate with a new key.
- If the credentials file is invalid or unreadable, run `cytetype setup --force` or remove it with `cytetype logout` before retrying setup.
- `cytetype logout` deletes only the local credentials. Revoke the API key from the dashboard if the server should stop accepting it.

See [CLI and Authentication](./cli.md) for the full command reference and credential locations.

## Annotation and Artifacts

- Rate limit responses include retry information; wait or provide your own LLM.
- Verify preprocessing: clustering and `rank_genes_groups` must be present.
- Make sure you have valid gene symbols in the AnnData object and are passing the correct gene symbols column name to parameter `gene_symbols_column`.
- If you are using a custom LLM, make sure you have the correct API key and base URL.
- For large datasets, load AnnData in backed mode (`sc.read_h5ad(..., backed="r")`) to reduce memory use during artifact generation.
- `run()` creates `vars.h5` and `obs.duckdb` before annotation. Use `cleanup_artifacts=True` if you do not want to keep these local files.
- If artifact building or uploading fails, `run()` will raise an error by default. Set `require_artifacts=False` to skip artifacts and continue with annotation only.
