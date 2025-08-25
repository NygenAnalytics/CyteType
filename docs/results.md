# Results

CyteType returns detailed annotations per cluster and stores them on your `AnnData` object.

## Access Results

### Using get_results() Method (Recommended)

The `get_results()` method is the recommended way to retrieve CyteType results. It provides intelligent result fetching with automatic fallback to API retrieval:

```python
# Helper
results = annotator.get_results()

import json
results = json.loads(adata.uns["cytetype_results"]["result"])
```

## Typical Fields (per cluster)
- annotation: primary cell type
- granularAnnotation: finer phenotype
- cellState: state/activation/malignancy
- ontologyTerm: Cell Ontology term (e.g., CL_0000127)
- confidence: model-assessed confidence
- supportingMarkers / conflictingMarkers
- missingExpression / unexpectedExpression
- corroboratingPapers: literature support

Example loop:
```python
for ann in results["annotations"]:
    print(ann["clusterId"], ann["annotation"], ann["ontologyTerm"], ann["confidence"]) 
```

For end-to-end examples and reports, see [Examples](./examples.md). 