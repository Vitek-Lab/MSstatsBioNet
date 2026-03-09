```r
#' Export network data with Cytoscape visualization
#' 
#' Convenience function that takes nodes and edges data directly and creates
#' both the configuration and HTML export in one step.
#' 
#' @inheritParams cytoscapeNetwork
#' @param filename Output HTML filename
#' @param ... Additional arguments passed to exportCytoscapeToHTML()
#' @export
#' @return Invisibly returns the file path of the created HTML file
exportNetworkToHTML <- function(nodes, edges, 
                                filename = "network_visualization.html",
                                displayLabelType = "id",
                                nodeFontSize = 12,
                                ...) {
    
    widget <- cytoscapeNetwork(nodes, edges, 
                               displayLabelType = displayLabelType,
                               nodeFontSize = nodeFontSize)
    
    htmlwidgets::saveWidget(
        widget,
        file = filename,
        selfcontained = TRUE
    )
    
    invisible(filename)
}

#' Preview network in browser
#' 
#' @description
#' Generates a temporary HTML file for the network visualization and opens it 
#' in the default web browser for quick preview.
#' 
#' @param nodes \code{data.frame} containing node information.
#' @param edges \code{data.frame} containing edge information.
#' @param displayLabelType \code{character} specifying the type of label to display on nodes. Default is \code{"id"}.
#' @param nodeFontSize \code{numeric} specifying the font size of node labels. Default is \code{12}.
#' 
#' @return Invisibly returns the file path of the temporary HTML file.
#' 
#' @examples
#' \dontrun{
#' nodes <- data.frame(id = c("A", "B", "C"))
#' edges <- data.frame(source = c("A", "B"), target = c("B", "C"))
#' previewNetworkInBrowser(nodes, edges)
#' }
#' 
#' @export
#' @importFrom utils browseURL
#' @inheritParams exportNetworkToHTML
previewNetworkInBrowser <- function(nodes, edges, 
                                    displayLabelType = "id",
                                    nodeFontSize = 12) {
    
    # Create temporary filename
    temp_file <- tempfile(fileext = ".html")
    
    # Export to temp file
    exportNetworkToHTML(nodes, edges, 
                        filename = temp_file,
                        displayLabelType = displayLabelType,
                        nodeFontSize = nodeFontSize)
    
    # Open in browser
    if (interactive()) {
        browseURL(temp_file)
        cat("Network opened in browser. Temporary file:", temp_file, "\n")
    }
    
    invisible(temp_file)
}
```