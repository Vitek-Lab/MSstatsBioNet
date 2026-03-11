# Preview network in browser

Generates a temporary HTML file for the network visualization and opens
it in the default web browser for quick preview.

## Usage

``` r
previewNetworkInBrowser(
  nodes,
  edges,
  displayLabelType = "id",
  nodeFontSize = 12
)
```

## Arguments

- nodes:

  Data frame with at minimum an `id` column. Optional columns: `logFC`
  (numeric), `hgncName` (character), `Site` (character,
  underscore-separated PTM site list).

- edges:

  Data frame with columns `source`, `target`, `interaction`. Optional:
  `site`, `evidenceLink`.

- displayLabelType:

  `"id"` (default) or `"hgncName"` – controls which column is used as
  the visible node label.

- nodeFontSize:

  Font size (px) for node labels. Default `12`.

## Value

Invisibly returns the file path of the temporary HTML file.

## Examples

``` r
if (FALSE) { # \dontrun{
nodes <- data.frame(id = c("A", "B", "C"))
edges <- data.frame(source = c("A", "B"), target = c("B", "C"))
previewNetworkInBrowser(nodes, edges)
} # }
```
