# Preview network in browser

Creates a temporary HTML file and opens it in the default web browser

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
