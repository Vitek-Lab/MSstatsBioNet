# MSstatsBioNet

[![Codecov test
coverage](https://codecov.io/github/Vitek-Lab/MSstatsBioNet/graph/badge.svg?token=SCPSPMTOEF)](https://codecov.io/github/Vitek-Lab/MSstatsBioNet)

This package provides a suite of functions to query various network
databases, filter queries & results, and visualize networks.

## Installation Instructions

To install this package on bioconductor, run the following command:

    if (!require("BiocManager", quietly = TRUE))
        install.packages("BiocManager")

    BiocManager::install("MSstatsBioNet")

You can install the development version of this package through Github:

    devtools::install_github("Vitek-Lab/MSstatsBioNet", build_vignettes = TRUE)

## Usage Examples

Here are some examples to help you get started with MSstatsBioNet:

### Annotate Protein Information

Use the `annotateProteinInfoFromIndra` function to annotate a data frame
with protein information from Indra.

``` r

library(MSstatsBioNet)

# Example data frame
df <- data.frame(Protein = c("CLH1_HUMAN"))

# Annotate protein information
annotated_df <- annotateProteinInfoFromIndra(df, "Uniprot_Mnemonic")
print(head(annotated_df))
```

### Visualize Networks with Cytoscape

Create an interactive network diagram using `cytoscapeNetwork`.

``` r

# Define nodes and edges
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

# Render the network
cytoscapeNetwork(nodes, edges)
```

### Export Network to HTML

Export your network visualization to an HTML file using
`exportNetworkToHTML`.

``` r

# Export the network to an HTML file
exportNetworkToHTML(nodes, edges, filename = "network.html")
```

### Retrieve Subnetwork from INDRA

Use `getSubnetworkFromIndra` to retrieve a subnetwork of protein
interactions from the INDRA database.

``` r

# Load example input data
input <- data.table::fread(system.file(
    "extdata/groupComparisonModel.csv",
    package = "MSstatsBioNet"
))

# Get subnetwork
subnetwork <- getSubnetworkFromIndra(input)
print(head(subnetwork$nodes))
print(head(subnetwork$edges))
```

### Preview Network in Browser

Quickly preview your network in a web browser using
`previewNetworkInBrowser`.

``` r

# Preview the network in a browser
previewNetworkInBrowser(nodes, edges)
```

### Integrate with Shiny

Use `cytoscapeNetworkOutput` and `renderCytoscapeNetwork` to integrate
network visualization into a Shiny app.

``` r

library(shiny)

ui <- fluidPage(
  cytoscapeNetworkOutput("cytoNetwork")
)

server <- function(input, output, session) {
  output$cytoNetwork <- renderCytoscapeNetwork({
    nodes <- data.frame(
      id = c("TP53", "MDM2", "CDKN1A"),
      logFC = c(1.5, -0.8, 2.1),
      stringsAsFactors = FALSE
    )
    edges <- data.frame(
      source = c("TP53", "MDM2"),
      target = c("MDM2", "TP53"),
      interaction = c("Activation", "Inhibition"),
      stringsAsFactors = FALSE
    )
    cytoscapeNetwork(nodes, edges)
  })
}

shinyApp(ui, server)
```

## License

This package is distributed under the
[Artistic-2.0](https://opensource.org/licenses/Artistic-2.0) license.
However, its dependencies may have different licenses.

Notably, INDRA is distributed under the [BSD
2-Clause](https://opensource.org/license/bsd-2-clause) license.
Furthermore, INDRA’s knowledge sources may have different licenses for
commercial applications. Please refer to the [INDRA
README](https://github.com/sorgerlab/indra?tab=readme-ov-file#indra-modules)
for more information on its knowledge sources and their associated
licenses.

## Databases Supported

- INDRA

## Filtering Options Supported

- P-Value Filter

## Visualization Options Supported

- Cytoscape Desktop
