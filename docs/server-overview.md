# Server Overview (High-Level)

CyteType client communicates with a hosted server that performs multi‑agent annotation of single‑cell clusters. The server is closed‑source; below is an operational overview to help you use the client effectively.

## Mental Model
- Submit once to start a job, then poll status and fetch results.
- Everything revolves around a `job_id`.
- You can re‑annotate a single cluster with feedback without re‑submitting the whole job.

## Key Endpoints
- POST `/annotate` — start a job. Optional `auth_token` for privacy
- GET `/status/{job_id}` — pending/processing/completed/failed + per‑cluster status
- GET `/results/{job_id}` — detailed results (summary + per‑cluster, with latest run)
- GET `/report/{job_id}` — HTML report shell backed by the same `/results`
- POST `/reannotate?job_id=...&cluster_id=...&feedback=...` — single‑cluster retry
- POST `/cluster_chat` — streaming Q&A on a cluster (SSE)

## Access Control
- If you submit with `auth_token`, the job is private; the same token is required for reads until you publish.
- Submissions without a token are public by default. Your job id is not indexed by search engines, so in effect it is private, unless you share the job id with someone.
- `GET /publish/{job_id}` to make a job public. This can be triggered via report.

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