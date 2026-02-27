#' Export network data with Cytoscape visualization
#' 
#' Convenience function that takes nodes and edges data directly and creates
#' both the configuration and HTML export in one step.
#' 
#' @param nodes Data frame with node information
#' @param edges Data frame with edge information  
#' @param filename Output HTML filename
#' @param displayLabelType Type of label to display ("id" or "hgncName")
#' @param nodeFontSize Font size for node labels (default: 12)
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
}

#' Preview network in browser
#' 
#' Creates a temporary HTML file and opens it in the default web browser
#' @export
#' @importFrom utils browseURL
#' @param nodes Data frame with node information
#' @param edges Data frame with edge information
#' @param displayLabelType Type of label to display ("id" or "hgncName")
#' @param ... Additional arguments passed to exportCytoscapeToHTML()
previewNetworkInBrowser <- function(nodes, edges, 
                                    displayLabelType = "id",
                                    nodeFontSize = 12,
                                    ...) {
    
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