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

- ...:

  Additional arguments passed to exportCytoscapeToHTML()

## Value

Invisibly returns the file path of the created HTML file
