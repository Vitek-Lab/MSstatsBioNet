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
  (numeric), `entityName` (character; may be semicolon-joined for
  multi-grounded rows), `entityId` (character), `Site` (character,
  underscore-separated PTM site list).

- edges:

  Data frame with columns `source`, `target`, `interaction`. Optional:
  `site`, `evidenceLink`.

- displayLabelType:

  `"id"` (default) or `"entityName"` – controls which column is used as
  the visible node label.

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
  logFC = c(1.5, -0.8, 2.1),
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
