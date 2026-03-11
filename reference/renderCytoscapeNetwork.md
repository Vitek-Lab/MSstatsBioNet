# Render a Cytoscape network in a Shiny application. This function is used to render a Cytoscape network visualization within a Shiny application.

Render a Cytoscape network in a Shiny application. This function is used
to render a Cytoscape network visualization within a Shiny application.

## Usage

``` r
renderCytoscapeNetwork(expr, env = parent.frame())
```

## Arguments

- expr:

  An expression that generates an HTML widget (or a
  [promise](https://rstudio.github.io/promises/) of an HTML widget).

- env:

  The environment in which to evaluate `expr`.

## Value

A rendered Cytoscape network widget for use in Shiny applications.

## Examples

``` r
if (FALSE) { # \dontrun{
library(shiny)
library(MSstatsBioNet)

ui <- fluidPage(
  cytoscapeNetworkOutput("cytoNetwork")
)

server <- function(input, output, session) {
  output$cytoNetwork <- renderCytoscapeNetwork({
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
  })
}

shinyApp(ui, server)
} # }
```
