# Server Overview (High-Level)

CyteType client communicates with a hosted server that performs multi‑agent annotation of single‑cell clusters. The server is closed‑source; below is an operational overview to help you use the client effectively.

## Mental Model
- Submit once to start a job, then poll status and fetch results.
- Everything revolves around a `job_id`.
- You can re‑annotate a single cluster with feedback without re‑submitting the whole job.

## Key Endpoints
- POST `/annotate`: start an authenticated job
- GET `/status/{job_id}`: pending/processing/completed/failed + per‑cluster status
- GET `/results/{job_id}`: detailed results (summary + per‑cluster, with latest run)
- GET `/report/{job_id}`: HTML report shell backed by the same `/results`
- POST `/reannotate?job_id=...&cluster_id=...&feedback=...`: single‑cluster retry
- POST `/cluster_chat`: streaming Q&A on a cluster (SSE)

## Authentication and Access
- The Python client requires a bearer token before uploading artifacts, submitting jobs, or fetching remote results.
- For local use, run `cytetype setup` and let the client load the saved credentials automatically.
- Saved credentials and direct tokens are tied to the selected API origin. Credentials for one server are not reused for another server.
- `cytetype view <job_id>` opens a report through the selected server's browser sign-in flow without putting the API key in the URL.

See [CLI and Authentication](./cli.md) for setup and credential handling.

## Rate Limits (typical defaults)
- Annotate: 5/day (Unlimited when a LLM is provided)
- Reannotate: 10/day
- Report/Status/Results: ~20/min
- Cluster chat: ~15/min

## What to Expect
- On submit, inputs and partial artifacts are saved; processing runs asynchronously per cluster.
- Results include annotations, ontology terms, evidence, literature, and usage summaries.
- History is preserved per cluster run; the server returns the latest by timestamp.

For more detail on how the multi‑agent system operates, see your deployment’s documentation or contact support.
