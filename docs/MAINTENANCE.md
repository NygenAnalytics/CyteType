# Documentation maintenance

Pending documentation work that is not tracked in the pages themselves. Keep HTML comments
out of `docs/source/`: the site publishes each page's source, page data and search index,
so comments become public.

## Update for new CyteTypeR releases

The pages describe CyteTypeR 0.6.0; list here the passages that a future CyteTypeR release changes.

## Screenshots to add

- `getting-started/how-cytetype-works.md`, section "What happens during a job": report live progress showing the six job stages and per-cluster status
- `getting-started/quickstart-r.md`, section "Plot the CyteType labels on the UMAP.": DimPlot of the PBMC 3k UMAP colored by cytetype_annotation_seurat_clusters
- `getting-started/sign-in.md`, section "Sign in with your browser": "Authorize CyteType CLI?" page with the account email and the Cancel and Authorize CLI buttons
- `getting-started/sign-in.md`, section "Remote servers, SSH and headless machines": dashboard API keys page with the "API key created" dialog and its Copy button
- `guides/prepare-data.md`, section "Cluster IDs and what you see in the report": report cluster card showing the User Label and its Match percentage
- `guides/running-jobs.md`, section "Follow a job": report while a job runs, with the progress stages and clusters still pending
- `guides/study-context.md`, section "Metadata tags": report Study Information panel with User-Provided Metadata, a linked tag and the client version
- `platform/chat.md`, section "Chat": Chat workspace with an attached run, a selected cluster, an answer with a figure and the Outputs panel open
- `platform/chat.md`, section "Outputs and downloads": Outputs panel with the output menu open (Add to chat, Download CSV, Provenance, Pin, Delete)
- `platform/dashboard.md`, section "Open the dashboard": dashboard Reports view with the sidebar, the workspace switcher and the My Reports table
- `platform/dashboard.md`, section "Share by email": Share Reports dialog with one active recipient and one Pending invite
- `platform/dashboard.md`, section "API keys": API keys page with the Create API key button, the Active and All filters and the key table
- `platform/dashboard.md`, section "Account": Account page with the Usage and limits card for the personal workspace
- `platform/mcp.md`, section "Claude (web, desktop and mobile)": Claude custom connector dialog with the CyteType name and server URL
- `platform/target-analysis.md`, section "Start an analysis": Target analysis dialog with the status line, the Analysis options (JSON) box and the Run analysis button
- `reference/cli.md`, section "`cytetype setup`": browser page "Authorize CyteType CLI?" with the Cancel and Authorize CLI buttons
- `results/report.md`, section "The CyteType report": full report with a cluster selected, showing the sidebar, the embedding with the cluster card and the right panel
- `results/report.md`, section "Follow a running job": report of a running job, with the Study Information progress tracker and pending clusters in the sidebar
- `results/report.md`, section "Left sidebar": left sidebar in Biological Order with the Filter menu open
- `results/report.md`, section "Gene search": embedding colored by a gene, with the gene legend and a hover value
- `results/report.md`, section "Cluster card": cluster card below the embedding, with the confidence badge, the Heterogeneous badge, the ontology link and the User Label match
- `results/report.md`, section "Toolbar": report toolbar with the Share report popover open
- `results/report.md`, section "Justification": Justification section with the annotation hierarchy and one expanded Selection Rounds entry
- `results/report.md`, section "Cell Type Evidence": Cell Type Evidence with the volcano plot and significant-gene table for one nearby cluster
- `results/report.md`, section "Quality Control": Quality Control with the cell classifier bars and one expanded reviewer in the Reviewer Panel
- `results/report.md`, section "Study summary": Dataset Overview with the Study Overview text and the Marker Gene Expression dot plot
