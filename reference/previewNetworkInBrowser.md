# Preview network in browser

Creates a temporary HTML file and opens it in the default web browser

## Usage

``` r
previewNetworkInBrowser(
  nodes,
  edges,
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

- displayLabelType:

  Type of label to display ("id" or "hgncName")

- ...:

  Additional arguments passed to exportCytoscapeToHTML()
