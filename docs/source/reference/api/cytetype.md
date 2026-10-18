---
myst:
  html_meta:
    description: Reference for the CyteType class, which prepares an AnnData object, submits an annotation job and stores the results in the AnnData object.
---

# CyteType class

`CyteType` annotates the clusters of one AnnData object. Creating it checks the inputs,
summarizes each cluster and writes the upload files `vars.h5` and `obs.duckdb`. `run()`
uploads the files, submits the job, waits for it to finish and adds the results to the
AnnData object. `get_results()` returns the stored results or fetches those of a finished
job, and `cleanup()` deletes the upload files.

```python
from cytetype import CyteType

annotator = CyteType(adata, group_key="leiden", rank_key="rank_genes_leiden")
adata = annotator.run(
    study_context="Peripheral blood mononuclear cells from a healthy human donor."
)
results = annotator.get_results()
annotator.cleanup()
```

Always use the AnnData object that `run()` returns. The
[Python quickstart](../getting-started/quickstart) runs these steps on a public dataset,
[Prepare your data](../guides/prepare-data) describes the input requirements, and
[Results in AnnData and Seurat](../results/results-in-your-data) describes every column
and key that `run()` adds.

```{eval-rst}
.. autoclass:: cytetype.CyteType
   :members: run, get_results, cleanup
```
