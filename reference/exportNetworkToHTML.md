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

  Data frame with node information

- edges:

  Data frame with edge information

- filename:

  Output HTML filename

- displayLabelType:

  Type of label to display ("id" or "hgncName")

- nodeFontSize:

  Font size for node labels (default: 12)

- ...:

  Additional arguments passed to exportCytoscapeToHTML()

## Value

Invisibly returns the file path of the created HTML file
