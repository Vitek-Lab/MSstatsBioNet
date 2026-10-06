# Render a Cytoscape network visualisation

Creates an interactive network diagram powered by Cytoscape.js and the
dagre layout algorithm. Nodes can carry log fold-change (logFC) values
which are mapped to a blue-grey-red colour gradient. PTM
(post-translational modification) site information is shown as small
satellite nodes and edge overlaps are surfaced as hover tooltips.

## Usage

``` r
cytoscapeNetwork(
  nodes,
  edges = data.frame(),
  displayLabelType = "id",
  nodeFontSize = 12,
  layoutOptions = NULL,
  width = NULL,
  height = NULL,
  elementId = NULL
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

- displayLabelType:

  `"id"` (default) or `"entity_name"`: which column is used as the
  visible node label. `"entityName"` is deprecated and is treated as
  `"entity_name"` with a warning.

- nodeFontSize:

  Font size (px) for node labels. Default `12`.

- layoutOptions:

  Named list of dagre layout options to override the defaults (e.g.
  `list(rankDir = "LR")`).

- width, height:

  Widget dimensions passed to
  [`createWidget`](https://rdrr.io/pkg/htmlwidgets/man/createWidget.html).

- elementId:

  Optional explicit HTML element id.

## Value

An `htmlwidget` object that renders in R Markdown, Shiny, or the RStudio
Viewer pane.

## Examples

``` r
if (FALSE) { # \dontrun{
nodes <- data.frame(
  id    = c("TP53", "MDM2", "CDKN1A"),
  logFC  = c(1.5, -0.8, 2.1),
  stringsAsFactors = FALSE
)
edges <- data.frame(
  source      = c("TP53",  "MDM2"),
  target      = c("MDM2",  "TP53"),
  interaction = c("Activation", "Inhibition"),
  stringsAsFactors = FALSE
)
cytoscapeNetwork(nodes, edges)
} # }
```
