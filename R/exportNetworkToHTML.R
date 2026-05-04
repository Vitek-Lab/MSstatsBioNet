```r
#' Export Network Data to HTML with Cytoscape Visualization
#' 
#' @description
#' This function exports network data consisting of nodes and edges to an HTML file
#' with a Cytoscape visualization. It combines the configuration and HTML export
#' into a single step, providing a convenient way to visualize networks.
#' 
#' @param nodes A \code{data.frame} containing the nodes of the network. Each row
#' represents a node with at least an \code{id} column.
#' @param edges A \code{data.frame} containing the edges of the network. Each row
#' represents an edge with at least \code{source} and \code{target} columns.
#' @param filename A \code{character} string specifying the output HTML filename.
#' Default is \code{"network_visualization.html"}.
#' @param displayLabelType A \code{character} string indicating the type of label
#' to display on nodes. Default is \code{"id"}.
#' @param nodeFontSize A \code{numeric} value specifying the font size of node labels.
#' Default is \code{12}.
#' @param ... Additional arguments passed to \code{exportCytoscapeToHTML()}.
#' 
#' @return Invisibly returns the file path of the created HTML file.
#' 
#' @examples
#' nodes <- data.frame(id = c("A", "B", "C"))
#' edges <- data.frame(source = c("A", "B"), target = c("B", "C"))
#' exportNetworkToHTML(nodes, edges, filename = "my_network.html")
#' 
#' @export
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