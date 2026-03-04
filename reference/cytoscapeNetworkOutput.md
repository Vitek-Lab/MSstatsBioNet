# Shiny output binding for cytoscapeNetwork

Creates a Shiny output binding for a Cytoscape network visualization,
allowing the network to be rendered within Shiny applications.

## Usage

``` r
cytoscapeNetworkOutput(outputId, width = "100%", height = "500px")
```

## Arguments

- outputId:

  output variable to read from

- width, height:

  Must be a valid CSS unit (like `"100%"`, `"400px"`, `"auto"`) or a
  number, which will be coerced to a string and have `"px"` appended.

## Value

A Shiny output binding for a Cytoscape network visualization.

## Examples

``` r
if (FALSE) { # \dontrun{
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
} # }
```
