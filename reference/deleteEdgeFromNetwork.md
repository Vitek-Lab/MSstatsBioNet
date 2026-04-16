# Delete an edge from a network edges data frame

Removes the row(s) from an edges data frame that match the given
`source`, `target`, and `interaction` values. This is the programmatic
counterpart of the interactive Ctrl+click / right-click edge deletion
available in
[`cytoscapeNetwork`](https://vitek-lab.github.io/MSstatsBioNet/reference/cytoscapeNetwork.md).

## Usage

``` r
deleteEdgeFromNetwork(edges, source, target, interaction)
```

## Arguments

- edges:

  Data frame with at minimum columns `source`, `target`, and
  `interaction`.

- source:

  Character. The source node identifier of the edge to remove.

- target:

  Character. The target node identifier of the edge to remove.

- interaction:

  Character. The interaction type of the edge to remove.

## Value

The `edges` data frame with the matching row(s) removed.

## Examples

``` r
edges <- data.frame(
  source      = c("TP53",  "MDM2",  "CDKN1A"),
  target      = c("MDM2",  "TP53",  "TP53"),
  interaction = c("Activation", "Inhibition", "Activation"),
  stringsAsFactors = FALSE
)
deleteEdgeFromNetwork(edges, "MDM2", "TP53", "Inhibition")
#>   source target interaction
#> 1   TP53   MDM2  Activation
#> 3 CDKN1A   TP53  Activation
```
