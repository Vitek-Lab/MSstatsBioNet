#' Render a Cytoscape network visualisation
#'
#' Creates an interactive network diagram powered by Cytoscape.js and the dagre
#' layout algorithm.  Nodes can carry log fold-change (logFC) values which are
#' mapped to a blue-grey-red colour gradient.  PTM (post-translational
#' modification) site information is shown as small satellite nodes and edge
#' overlaps are surfaced as hover tooltips.
#'
#' @param nodes       Data frame with at minimum an \code{id} column.  Optional
#'                    columns: \code{logFC} (numeric), \code{hgncName}
#'                    (character), \code{Site} (character, underscore-separated
#'                    PTM site list).
#' @param edges       Data frame with columns \code{source}, \code{target},
#'                    \code{interaction}.  Optional: \code{site},
#'                    \code{evidenceLink}.
#' @param displayLabelType \code{"id"} (default) or \code{"hgncName"} –
#'                    controls which column is used as the visible node label.
#' @param nodeFontSize Font size (px) for node labels.  Default \code{12}.
#' @param layoutOptions Named list of dagre layout options to override the
#'                    defaults (e.g. \code{list(rankDir = "LR")}).
#' @param networkId    Optional string used as the localStorage key for saving
#'                    the network state (node positions and deleted edges)
#'                    between page loads.  If \code{NULL} (the default), a
#'                    stable identifier is derived from the node and edge data
#'                    so the correct state is restored when the same network is
#'                    re-opened.  Supply an explicit value when you export
#'                    multiple different networks to the same file and want to
#'                    keep their saved states separate.
#' @param width,height Widget dimensions passed to
#'                    \code{\link[htmlwidgets]{createWidget}}.
#' @param elementId   Optional explicit HTML element id.
#'
#' @return An \code{htmlwidget} object that renders in R Markdown, Shiny, or
#'   the RStudio Viewer pane.
#'
#' @examples
#' \dontrun{
#' nodes <- data.frame(
#'   id    = c("TP53", "MDM2", "CDKN1A"),
#'   logFC = c(1.5, -0.8, 2.1),
#'   stringsAsFactors = FALSE
#' )
#' edges <- data.frame(
#'   source      = c("TP53",  "MDM2"),
#'   target      = c("MDM2",  "TP53"),
#'   interaction = c("Activation", "Inhibition"),
#'   stringsAsFactors = FALSE
#' )
#' cytoscapeNetwork(nodes, edges)
#' }
#'
#' @importFrom htmlwidgets createWidget
#' @importFrom grDevices colorRamp rgb
#' @export
cytoscapeNetwork <- function(nodes,
                             edges            = data.frame(),
                             displayLabelType = "id",
                             nodeFontSize     = 12,
                             layoutOptions    = NULL,
                             networkId        = NULL,
                             width            = NULL,
                             height           = NULL,
                             elementId        = NULL) {
    
    # Validate inputs
    if (!is.data.frame(nodes) || !("id" %in% names(nodes))) {
        stop("`nodes` must be a data frame with at least an `id` column.")
    }
    if (!is.data.frame(edges)) {
        stop("`edges` must be a data frame.")
    }
    required_edge_cols <- c("source", "target", "interaction")
    if (nrow(edges) > 0 && !all(required_edge_cols %in% names(edges))) {
        stop("`edges` must contain columns: source, target, interaction.")
    }
    
    # Build layout config
    default_layout <- list(
        name          = "dagre",
        rankDir       = "TB",
        animate       = TRUE,
        fit           = TRUE,
        padding       = 30,
        spacingFactor = 1.5,
        nodeSep       = 50,
        edgeSep       = 20,
        rankSep       = 80
    )
    layout <- default_layout
    if (!is.null(layoutOptions)) {
        for (nm in names(layoutOptions)) layout[[nm]] <- layoutOptions[[nm]]
    }
    
    # Build element list
    elements <- .buildElements(nodes, edges, displayLabelType)

    # Derive a stable networkId from the data if none was supplied.
    # This ensures the correct localStorage state is restored when the same
    # network is re-opened, while different networks get different keys.
    if (is.null(networkId)) {
        node_str <- paste(sort(as.character(nodes$id)), collapse = "|")
        edge_str <- if (nrow(edges) > 0) {
            paste(sort(paste(edges$source, edges$target, edges$interaction,
                             sep = "-")),
                  collapse = "|")
        } else ""
        combined  <- paste0(node_str, "__", edge_str)
        ints      <- utf8ToInt(combined)
        hash_val  <- Reduce(function(acc, x) (acc * 31L + x) %% 16777216L,
                            ints, 0L)
        networkId <- sprintf("net%06x", hash_val)
    }

    # Package everything for the JS side
    x <- list(
        elements       = elements,
        layout         = layout,
        node_font_size = nodeFontSize,
        network_id     = networkId
    )
    
    htmlwidgets::createWidget(
        name      = "cytoscapeNetwork",
        x         = x,
        width     = width,
        height    = height,
        package   = "MSstatsBioNet",
        elementId = elementId
    )
}


# ── Shiny helpers ───────────────────────────────────────────────────────────

#' Shiny output binding for cytoscapeNetwork
#'
#' Creates a Shiny output binding for a Cytoscape network visualization, allowing
#' the network to be rendered within Shiny applications.
#'
#' @return A Shiny output binding for a Cytoscape network visualization.
#'
#' @examples
#' \dontrun{
#' library(shiny)
#' 
#' ui <- fluidPage(
#'   cytoscapeNetworkOutput("cytoNetwork")
#' )
#' 
#' server <- function(input, output, session) {
#'   output$cytoNetwork <- renderCytoscapeNetwork({
#'     nodes <- data.frame(
#'       id = c("TP53", "MDM2", "CDKN1A"),
#'       logFC = c(1.5, -0.8, 2.1),
#'       stringsAsFactors = FALSE
#'     )
#'     edges <- data.frame(
#'       source = c("TP53", "MDM2"),
#'       target = c("MDM2", "TP53"),
#'       interaction = c("Activation", "Inhibition"),
#'       stringsAsFactors = FALSE
#'     )
#'     cytoscapeNetwork(nodes, edges)
#'   })
#' }
#' 
#' shinyApp(ui, server)
#' }
#'
#' @importFrom htmlwidgets shinyWidgetOutput
#' @inheritParams htmlwidgets::shinyWidgetOutput
#' @export
cytoscapeNetworkOutput <- function(outputId,
                                   width  = "100%",
                                   height = "500px") {
    htmlwidgets::shinyWidgetOutput(
        outputId = outputId,
        name     = "cytoscapeNetwork",
        width    = width,
        height   = height,
        package  = "MSstatsBioNet"
    )
}

#' Render a Cytoscape network in a Shiny application.
#' This function is used to render a Cytoscape network visualization within a Shiny application.
#' 
#' @importFrom htmlwidgets shinyRenderWidget createWidget
#' @inheritParams htmlwidgets::shinyRenderWidget
#'
#' @return A rendered Cytoscape network widget for use in Shiny applications.
#'
#' @examples
#' \dontrun{
#' library(shiny)
#' library(MSstatsBioNet)
#'
#' ui <- fluidPage(
#'   cytoscapeNetworkOutput("cytoNetwork")
#' )
#'
#' server <- function(input, output, session) {
#'   output$cytoNetwork <- renderCytoscapeNetwork({
#'     nodes <- data.frame(
#'       id    = c("TP53", "MDM2", "CDKN1A"),
#'       logFC = c(1.5, -0.8, 2.1),
#'       stringsAsFactors = FALSE
#'     )
#'     edges <- data.frame(
#'       source      = c("TP53",  "MDM2"),
#'       target      = c("MDM2",  "TP53"),
#'       interaction = c("Activation", "Inhibition"),
#'       stringsAsFactors = FALSE
#'     )
#'     cytoscapeNetwork(nodes, edges)
#'   })
#' }
#'
#' shinyApp(ui, server)
#' }
#' @export
renderCytoscapeNetwork <- function(expr, env = parent.frame()) {
    expr <- substitute(expr)
    htmlwidgets::shinyRenderWidget(
        expr    = expr,
        outputFunction = cytoscapeNetworkOutput,
        env     = env,
        quoted  = TRUE
    )
}