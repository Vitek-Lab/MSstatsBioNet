# Create visualization of network

Use results from INDRA to generate a visualization of the a network on
Cytoscape Desktop. Note that the Cytoscape Desktop app must be open for
this function to work.

## Usage

``` r
visualizeNetworks(
  nodes,
  edges,
  pvalueCutoff = 0.05,
  logfcCutoff = 0.5,
  node_label_column = "id",
  main_targets = c()
)
```

## Arguments

- nodes:

  dataframe of nodes consisting of columns id (chararacter), pvalue
  (number), logFC (number)

- edges:

  dataframe of edges consisting of columns source (character), target
  (character), interaction (character), evidenceCount (number),
  evidenceLink (character)

- pvalueCutoff:

  p-value cutoff for coloring significant proteins. Default is 0.05

- logfcCutoff:

  log fold change cutoff for coloring significant proteins. Default is
  0.5

- node_label_column:

  The column of the nodes dataframe to use as the node label. Default is
  "id". "hgncName" can be used for gene name.

- main_targets:

  character vector of main targets to stand-out with a different node
  shape. Default is an empty vector c(). IDs of main targets should
  match the column used by the node_label_column parameter.

## Value

cytoscape visualization of subnetwork

## Examples

``` r
input <- data.table::fread(system.file(
    "extdata/groupComparisonModel.csv",
    package = "MSstatsBioNet"
))
subnetwork <- getSubnetworkFromIndra(input)
#> Warning: NOTICE: This function includes third-party software components
#>         that are licensed under the BSD 2-Clause License. Please ensure to
#>         include the third-party licensing agreements if redistributing this
#>         package or utilizing the results based on this package.
#>         See the LICENSE file for more details.
visualizeNetworks(subnetwork$nodes, subnetwork$edges)
#> Warning: Visualization is not available in non-interactive mode.
```
