# R/cytoscapeNetwork.R
#
# htmlwidgets binding for the Cytoscape network visualisation.
# The heavy-lifting JS lives in inst/htmlwidgets/cytoscapeNetwork.js.
# This file is responsible for:
#   1. Pre-processing nodes / edges in R (colour mapping, element serialisation)
#   2. Calling htmlwidgets::createWidget() to hand everything to the JS side.

# ── Internal helpers (not exported) ────────────────────────────────────────

#' Map logFC values to a blue-grey-red colour palette
#' @keywords internal
#' @noRd
.mapLogFCToColor <- function(logFC_values) {
  colors <- c("#ADD8E6", "#ADD8E6", "#D3D3D3", "#FFA590", "#FFA590")

  if (all(is.na(logFC_values)) ||
      length(unique(logFC_values[!is.na(logFC_values)])) <= 1) {
    return(rep("#D3D3D3", length(logFC_values)))
  }

  default_max     <- 2
  max_logFC       <- max(c(abs(logFC_values), default_max), na.rm = TRUE)
  min_logFC       <- -max_logFC
  color_map       <- grDevices::colorRamp(colors)
  normalized      <- (logFC_values - min_logFC) / (max_logFC - min_logFC)
  normalized[is.na(normalized)] <- 0.5
  rgb_colors      <- color_map(normalized)
  grDevices::rgb(rgb_colors[, 1], rgb_colors[, 2], rgb_colors[, 3],
                 maxColorValue = 255)
}

#' Safely escape a string for embedding in a JS single-quoted literal
#' @keywords internal
#' @noRd
.escJS <- function(x) {
  if (is.null(x)) return("")
  x <- as.character(x)
  x <- gsub("\\\\", "\\\\\\\\", x)
  x <- gsub("'",    "\\\\'",    x)
  x <- gsub("\r",   "\\\\r",    x)
  x <- gsub("\n",   "\\\\n",    x)
  x
}

#' Relationship properties lookup
#' @keywords internal
#' @noRd
.relProps <- function() {
  list(
    complex = list(
      types       = "Complex",
      color       = "#8B4513",
      style       = "solid",
      arrow       = "none",
      width       = 4,
      consolidate = "undirected"
    ),
    regulatory = list(
      types       = c("Inhibition", "Activation", "IncreaseAmount", "DecreaseAmount"),
      colors      = list(Inhibition     = "#FF4444",
                         Activation     = "#44AA44",
                         IncreaseAmount = "#4488FF",
                         DecreaseAmount = "#FF8844"),
      style       = "solid",
      arrow       = "triangle",
      width       = 3,
      consolidate = "bidirectional"
    ),
    phosphorylation = list(
      types       = "Phosphorylation",
      color       = "#9932CC",
      style       = "dashed",
      arrow       = "triangle",
      width       = 2,
      consolidate = "directed"
    ),
    other = list(
      color       = "#666666",
      style       = "dotted",
      arrow       = "triangle",
      width       = 2,
      consolidate = "directed"
    )
  )
}

#' Classify an interaction string into a relationship category
#' @keywords internal
#' @noRd
.classify <- function(interaction) {
  props <- .relProps()
  for (cat_name in names(props)) {
    if (!is.null(props[[cat_name]]$types) &&
        interaction %in% props[[cat_name]]$types) {
      return(cat_name)
    }
  }
  "other"
}

#' Retrieve edge colour / style / arrow / width
#' @keywords internal
#' @noRd
.edgeStyle <- function(interaction, category, edge_type) {
  props <- .relProps()
  p     <- if (category %in% names(props)) props[[category]] else props$other

  color <- if (category == "regulatory" && !is.null(p$colors)) {
    base <- sub(" \\(bidirectional\\)", "", interaction)
    if (base %in% names(p$colors)) p$colors[[base]] else "#666666"
  } else {
    p$color
  }

  arrow <- switch(edge_type,
    undirected   = "none",
    bidirectional = "triangle",
    p$arrow
  )

  list(color = color, style = p$style, arrow = arrow, width = p$width)
}

#' Aggregate PTM overlap between edge targets and node Site columns
#' @keywords internal
#' @noRd
.ptmOverlap <- function(edges, nodes) {
  if (nrow(edges) == 0 || is.null(nodes)) return(setNames(character(0), character(0)))

  edges$edge_key <- paste(edges$source, edges$target, edges$interaction, sep = "-")
  unique_keys    <- unique(edges$edge_key)
  result         <- setNames(character(length(unique_keys)), unique_keys)

  for (key in unique_keys) {
    sub_edges <- edges[edges$edge_key == key, ]
    all_sites <- c()

    for (i in seq_len(nrow(sub_edges))) {
      e <- sub_edges[i, ]
      if (!is.na(e$target) && "site" %in% names(e) && !is.na(e$site)) {
        tnodes <- nodes[nodes$id == e$target, ]
        if (nrow(tnodes) > 0 && "Site" %in% names(tnodes)) {
          edge_sites <- trimws(unlist(strsplit(as.character(e$site), "[,;|]")))
          for (j in seq_len(nrow(tnodes))) {
            if (!is.na(tnodes$Site[j])) {
              node_sites   <- trimws(unlist(strsplit(as.character(tnodes$Site[j]), "_")))
              overlap      <- intersect(edge_sites, node_sites)
              overlap      <- overlap[overlap != "" & !is.na(overlap)]
              all_sites    <- c(all_sites, overlap)
            }
          }
        }
      }
    }

    u <- unique(all_sites[all_sites != "" & !is.na(all_sites)])
    result[key] <- if (length(u) == 0) {
      ""
    } else if (length(u) == 1) {
      paste0("Overlapping PTM site: ",  u)
    } else {
      paste0("Overlapping PTM sites: ", paste(u, collapse = ", "))
    }
  }
  result
}

#' Consolidate bidirectional / undirected edges
#' @keywords internal
#' @noRd
.consolidateEdges <- function(edges, nodes = NULL) {
  if (nrow(edges) == 0) return(edges)

  ptm_map      <- .ptmOverlap(edges, nodes)
  props        <- .relProps()
  consolidated <- list()
  processed    <- c()

  for (i in seq_len(nrow(edges))) {
    e        <- edges[i, ]
    pair_key <- paste(sort(c(e$source, e$target)), e$interaction, collapse = "-")
    if (pair_key %in% processed) next

    cat        <- .classify(e$interaction)
    rev_edges  <- edges[edges$source == e$target &
                        edges$target == e$source &
                        edges$interaction == e$interaction, ]
    con_type   <- props[[cat]]$consolidate
    edge_key   <- paste(e$source, e$target, e$interaction, sep = "-")
    ptm_txt    <- if (edge_key %in% names(ptm_map)) ptm_map[[edge_key]] else ""

    if (nrow(rev_edges) > 0 && con_type %in% c("undirected", "bidirectional")) {
      new_interaction <- if (con_type == "undirected") e$interaction else
        paste(e$interaction, "(bidirectional)")
      new_edge <- data.frame(source      = e$source,
                             target      = e$target,
                             interaction = new_interaction,
                             edge_type   = if (con_type == "undirected") "undirected" else "bidirectional",
                             category    = cat,
                             ptm_overlap = ptm_txt,
                             stringsAsFactors = FALSE)
      for (col in setdiff(names(e), c("source", "target", "interaction"))) {
        new_edge[[col]] <- e[[col]]
      }
      key <- paste(e$source, e$target, new_interaction, sep = "-")
      consolidated[[key]] <- new_edge
      processed <- c(processed, pair_key)
    } else {
      de            <- e
      de$edge_type  <- "directed"
      de$category   <- cat
      de$ptm_overlap <- ptm_txt
      key           <- paste(e$source, e$target, e$interaction, sep = "-")
      consolidated[[key]] <- de
    }
  }

  if (length(consolidated) > 0) {
    result           <- do.call(rbind, consolidated)
    rownames(result) <- NULL
    result
  } else {
    edges[0, ]
  }
}

#' Build the list of Cytoscape element objects (nodes + edges)
#'
#' Returns a list of named lists — jsonlite will serialise them cleanly.
#' @keywords internal
#' @noRd
.buildElements <- function(nodes, edges, display_label_type = "id") {
  # ── node colours ──────────────────────────────────────────────────────
  node_colors <- if ("logFC" %in% names(nodes)) {
    .mapLogFCToColor(nodes$logFC)
  } else {
    rep("#D3D3D3", nrow(nodes))
  }

  label_col <- if (display_label_type == "hgncName" &&
                   "hgncName" %in% names(nodes)) "hgncName" else "id"

  has_ptm_sites <- if ("Site" %in% names(nodes)) {
    unique(nodes$id[!is.na(nodes$Site) & trimws(nodes$Site) != ""])
  } else {
    character(0)
  }

  elements        <- list()
  emitted_prots   <- character(0)
  emitted_cpds    <- character(0)
  emitted_ptm_n   <- character(0)
  emitted_ptm_e   <- character(0)

  for (i in seq_len(nrow(nodes))) {
    row       <- nodes[i, , drop = FALSE]
    color     <- node_colors[i]
    has_site  <- "Site" %in% names(nodes) &&
                 !is.na(row$Site) && trimws(row$Site) != ""

    display_label <- if (label_col == "hgncName" &&
                         !is.na(row$hgncName) && row$hgncName != "")
      row$hgncName else row$id

    needs_compound <- row$id %in% has_ptm_sites
    compound_id    <- paste0(row$id, "__compound__")

    # Compound container
    if (needs_compound && !(compound_id %in% emitted_cpds)) {
      elements <- c(elements, list(
        list(data = list(id        = compound_id,
                         node_type = "compound"))
      ))
      emitted_cpds <- c(emitted_cpds, compound_id)
    }

    # Protein node
    if (!(row$id %in% emitted_prots)) {
      nd <- list(id        = row$id,
                 label     = display_label,
                 color     = color,
                 node_type = "protein")
      if (needs_compound) nd$parent <- compound_id
      elements <- c(elements, list(list(data = nd)))
      emitted_prots <- c(emitted_prots, row$id)
    }

    # PTM child nodes + attachment edges
    if (has_site) {
      sites <- unique(trimws(unlist(strsplit(as.character(row$Site), "[_,;|]"))))
      sites <- sites[sites != ""]

      for (site in sites) {
        ptm_nid <- paste0(row$id, "__ptm__", site)
        if (!(ptm_nid %in% emitted_ptm_n)) {
          elements <- c(elements, list(list(data = list(
            id             = ptm_nid,
            label          = site,
            color          = color,
            parent_protein = row$id,
            parent         = compound_id,
            node_type      = "ptm"
          ))))
          emitted_ptm_n <- c(emitted_ptm_n, ptm_nid)
        }

        ptm_eid <- paste0(row$id, "__ptm_edge__", site)
        if (!(ptm_eid %in% emitted_ptm_e)) {
          elements <- c(elements, list(list(data = list(
            id          = ptm_eid,
            source      = row$id,
            target      = ptm_nid,
            edge_type   = "ptm_attachment",
            category    = "ptm_attachment",
            interaction = "",
            color       = color,
            line_style  = "dotted",
            arrow_shape = "none",
            width       = 1.5,
            tooltip     = ""
          ))))
          emitted_ptm_e <- c(emitted_ptm_e, ptm_eid)
        }
      }
    }
  }

  # ── edges ─────────────────────────────────────────────────────────────
  if (!is.null(edges) && nrow(edges) > 0) {
    con <- .consolidateEdges(edges, nodes)

    for (i in seq_len(nrow(con))) {
      row  <- con[i, ]
      sty  <- .edgeStyle(row$interaction, row$category, row$edge_type)
      eid  <- paste(row$source, row$target, row$interaction, sep = "-")
      elink <- if ("evidenceLink" %in% names(row)) {
        ev <- row$evidenceLink
        if (is.na(ev) || ev == "NA") "" else as.character(ev)
      } else ""

      elements <- c(elements, list(list(data = list(
        id          = eid,
        source      = row$source,
        target      = row$target,
        interaction = row$interaction,
        edge_type   = row$edge_type,
        category    = row$category,
        evidenceLink = elink,
        color       = sty$color,
        line_style  = sty$style,
        arrow_shape = sty$arrow,
        width       = sty$width,
        tooltip     = if (!is.null(row$ptm_overlap)) row$ptm_overlap else ""
      ))))
    }
  }

  elements
}


# ── Public API ──────────────────────────────────────────────────────────────

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
                             edges         = data.frame(),
                             displayLabelType = "id",
                             nodeFontSize  = 12,
                             layoutOptions = NULL,
                             width         = NULL,
                             height        = NULL,
                             elementId     = NULL) {

  # Validate inputs
  if (!is.data.frame(nodes) || !("id" %in% names(nodes))) {
    stop("`nodes` must be a data frame with at least an `id` column.")
  }
  if (!is.data.frame(edges)) edges <- data.frame()

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

  # Package everything for the JS side
  x <- list(
    elements       = elements,
    layout         = layout,
    node_font_size = nodeFontSize
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
    package  = "cytoscapeNetwork"
  )
}

#' Shiny render binding for cytoscapeNetwork
#' @inheritParams htmlwidgets::shinyRenderWidget
#' @export
renderCytoscapeNetwork <- function(expr, env = parent.frame(), quoted = FALSE) {
  if (!quoted) expr <- substitute(expr)
  htmlwidgets::shinyRenderWidget(
    expr    = expr,
    outputFunction = cytoscapeNetworkOutput,
    env     = env,
    quoted  = TRUE
  )
}
