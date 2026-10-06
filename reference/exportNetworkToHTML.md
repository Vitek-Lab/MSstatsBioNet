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

  Data frame with at minimum an `id` column. Optional columns: `logFC`
  (numeric), `entity_name` (character; may be semicolon-joined for
  multi-grounded rows), `entity_id` (character), `site` (character,
  underscore-separated PTM site list), and the node status columns of
  [`validate_network()`](https://vitek-lab.github.io/MSstatsBioNet/reference/validate_network.md):
  `entity_type` sets the node shape (e.g. a hexagon for `"metabolite"`);
  `measured = FALSE` draws a node with no fill and a dashed grey border
  ("not in input data"); `included_in_query = FALSE` fades it. A node
  with an `NA` `logFC` has no fill. A protein with `site` rows takes its
  colour from its row with no `site`; without one, only its sites are
  coloured. Missing status columns count as `TRUE`.

- edges:

  Data frame with columns `source`, `target`, `interaction`. Optional:
  `directed` (an edge with `FALSE` has no arrow and is drawn once for
  both directions; without the column, `Complex` and `Association` edges
  are undirected), `site`, `evidence_url`, opened when an edge is
  clicked.

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
