# Export network data with Cytoscape visualization

Convenience function that takes nodes and edges data directly and
creates both the configuration and HTML export in one step.

## Usage

``` r
exportNetworkToHTML(
  nodes,
  edges,
  filename = "network_visualization.html",
  displayLabelType = "id",
  nodeFontSize = 12,
  ...
)
```

## Arguments

- nodes:

  Data frame with at minimum an `id` column. Optional columns: `log2FC`
  (numeric), `entity_name` (character; may be semicolon-joined for
  multi-grounded rows), `entity_id` (character), `site` (character,
  underscore-separated PTM site list).

- edges:

  Data frame with columns `source`, `target`, `interaction`. Optional:
  `site`, `evidence_url`, opened when an edge is clicked.

- filename:

  Output HTML filename

- displayLabelType:

  `"id"` (default) or `"entity_name"`: which column is used as the
  visible node label. `"entityName"` is deprecated and is treated as
  `"entity_name"` with a warning.

- nodeFontSize:

  Font size (px) for node labels. Default `12`.

- ...:

  Additional arguments passed to exportCytoscapeToHTML()

## Value

Invisibly returns the file path of the created HTML file
