#' Node fill for a logFC of zero, and for nodes when there is no logFC column
#' @keywords internal
#' @noRd
NEUTRAL_NODE_COLOR <- "#D3D3D3"

#' Node fill for an NA logFC: distinct from the zero-fold-change grey
#' @keywords internal
#' @noRd
NO_DATA_NODE_COLOR <- "#FFFFFF"

#' Node shape for each nodes$entity_type
#'
#' Types that are not listed, and NA, get the protein shape. A PTM protein's
#' rows can be ptm_site rows; the protein node still gets the protein shape.
#' @keywords internal
#' @noRd
NODE_SHAPES <- c(protein    = "round-rectangle",
                 gene       = "round-rectangle",
                 transcript = "round-rectangle",
                 ptm_site   = "round-rectangle",
                 metabolite = "hexagon",
                 lipid      = "hexagon",
                 drug       = "diamond",
                 complex    = "octagon",
                 family     = "barrel",
                 other      = "round-rectangle")

#' Map logFC values to a blue-grey-red colour palette
#'
#' \code{NA} maps to \code{NO_DATA_NODE_COLOR}, not to the grey of a zero
#' fold change.
#' @importFrom grDevices colorRamp rgb
#' @keywords internal
#' @noRd
.mapLogFCToColor <- function(logFC_values) {
    colors <- c("#ADD8E6", "#ADD8E6", NEUTRAL_NODE_COLOR, "#FFA590", "#FFA590")
    is_missing <- is.na(logFC_values)
    result <- rep(NO_DATA_NODE_COLOR, length(logFC_values))

    if (length(unique(logFC_values[!is_missing])) <= 1) {
        result[!is_missing] <- NEUTRAL_NODE_COLOR
        return(result)
    }
    
    is_pos_inf <- is.infinite(logFC_values) & logFC_values > 0
    is_neg_inf <- is.infinite(logFC_values) & logFC_values < 0
    
    finite_values   <- logFC_values[is.finite(logFC_values)]
    default_max     <- 2
    max_logFC       <- max(c(abs(finite_values), default_max), na.rm = TRUE)
    min_logFC       <- -max_logFC
    
    logFC_values[is_pos_inf] <-  max_logFC
    logFC_values[is_neg_inf] <-  min_logFC
    
    color_map       <- grDevices::colorRamp(colors)
    normalized      <- (logFC_values[!is_missing] - min_logFC) /
        (max_logFC - min_logFC)
    rgb_colors      <- color_map(normalized)
    result[!is_missing] <- grDevices::rgb(rgb_colors[, 1], rgb_colors[, 2],
                                          rgb_colors[, 3], maxColorValue = 255)
    result
}

#' How a node is drawn, from the node status columns
#'
#' \describe{
#'   \item{latent}{\code{measured = FALSE}: not in the input data}
#'   \item{no_logfc}{in the input, but \code{logFC} is \code{NA}}
#'   \item{not_queried}{in the input, but \code{included_in_query = FALSE}}
#'   \item{measured}{everything else}
#' }
#' A missing \code{measured} or \code{included_in_query} column counts as
#' \code{TRUE}, and a missing \code{logFC} column gives no
#' \code{"no_logfc"} rows, so nodes built by hand render as before.
#' \code{.buildElements()} adds \code{"sites_only"} for a PTM protein that
#' has no protein-level row.
#' @param nodes nodes data.frame
#' @return character vector, one status per row
#' @keywords internal
#' @noRd
.node_display_status <- function(nodes) {
    measured <- .get_column_or_default(nodes, "measured", TRUE) %in%
        c(TRUE, NA)
    in_query <- .get_column_or_default(nodes, "included_in_query", TRUE) %in%
        c(TRUE, NA)
    no_logfc <- if ("logFC" %in% names(nodes)) {
        is.na(nodes$logFC)
    } else {
        rep(FALSE, nrow(nodes))
    }
    ifelse(!measured, "latent",
           ifelse(no_logfc, "no_logfc",
                  ifelse(!in_query, "not_queried", "measured")))
}

#' Node shape for each entity type
#' @param entity_types character vector of nodes$entity_type values
#' @return character vector of Cytoscape shape names
#' @keywords internal
#' @noRd
.node_shape <- function(entity_types) {
    shapes <- unname(NODE_SHAPES[as.character(entity_types)])
    shapes[is.na(shapes)] <- NODE_SHAPES[["protein"]]
    shapes
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
#'
#' Each statement type in \code{INTERACTION_TYPES} belongs to one category.
#' Whether an edge has an arrow comes from its \code{directed} column, not
#' from the category.
#' @keywords internal
#' @noRd
.relProps <- function() {
    list(
        complex = list(
            types       = c("Complex", "Association"),
            colors      = list(Complex     = "#8B4513",
                               Association = "#C08A5B"),
            style       = "solid",
            arrow       = "none",
            width       = 4
        ),
        regulatory = list(
            types       = c("Inhibition", "Activation", "IncreaseAmount",
                            "DecreaseAmount", "Regulation", "Influence",
                            "Gef", "Gap", "GtpActivation"),
            colors      = list(Inhibition     = "#FF4444",
                               Activation     = "#44AA44",
                               IncreaseAmount = "#4488FF",
                               DecreaseAmount = "#FF8844"),
            style       = "solid",
            arrow       = "triangle",
            width       = 3
        ),
        phosphorylation = list(
            types       = c("Phosphorylation", "Dephosphorylation",
                            "Autophosphorylation", "Transphosphorylation"),
            color       = "#9932CC",
            style       = "dashed",
            arrow       = "triangle",
            width       = 2
        ),
        modification = list(
            types       = c("Modification",
                            "Ubiquitination", "Deubiquitination",
                            "Sumoylation", "Desumoylation",
                            "Hydroxylation", "Dehydroxylation",
                            "Acetylation", "Deacetylation",
                            "Glycosylation", "Deglycosylation",
                            "Farnesylation", "Defarnesylation",
                            "Geranylgeranylation", "Degeranylgeranylation",
                            "Palmitoylation", "Depalmitoylation",
                            "Myristoylation", "Demyristoylation",
                            "Ribosylation", "Deribosylation",
                            "Methylation", "Demethylation"),
            color       = "#20898B",
            style       = "dashed",
            arrow       = "triangle",
            width       = 2
        ),
        other = list(
            types       = c("Conversion", "Translocation"),
            color       = "#666666",
            style       = "dotted",
            arrow       = "triangle",
            width       = 2
        )
    )
}

#' Classify an interaction string into a relationship category
#' @keywords internal
#' @noRd
.classify <- function(interaction) {
    if (is.null(interaction) || is.na(interaction) || !nzchar(trimws(as.character(interaction)))) {
        return("other")
    }
    interaction <- as.character(interaction)
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
    
    color <- if (!is.null(p$colors)) {
        if (interaction %in% names(p$colors)) p$colors[[interaction]] else "#666666"
    } else {
        p$color
    }

    arrow <- switch(edge_type,
                    undirected = "none",
                    p$arrow
    )
    
    list(color = color, style = p$style, arrow = arrow, width = p$width)
}

#' Aggregate PTM overlap between edge targets and node site columns
#' @keywords internal
#' @importFrom stats setNames
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
                if (nrow(tnodes) > 0 && "site" %in% names(tnodes)) {
                    edge_sites <- trimws(unlist(strsplit(as.character(e$site), "[,;|]")))
                    for (j in seq_len(nrow(tnodes))) {
                        if (!is.na(tnodes$site[j])) {
                            node_sites   <- trimws(unlist(strsplit(as.character(tnodes$site[j]), "_")))
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

#' Whether each edge is directed
#'
#' Reads the contract's \code{directed} column. Edges without it, or with
#' \code{NA}, are undirected when their statement type is in
#' \code{UNDIRECTED_INTERACTION_TYPES}.
#' @param edges edges data.frame
#' @return logical vector, one per edge
#' @keywords internal
#' @noRd
.edge_is_directed <- function(edges) {
    by_type <- !edges$interaction %in% UNDIRECTED_INTERACTION_TYPES
    if (!"directed" %in% names(edges)) {
        return(by_type)
    }
    directed <- as.logical(edges$directed)
    ifelse(is.na(directed), by_type, directed)
}

#' Consolidate undirected edges
#'
#' An undirected edge is drawn without an arrow, and its reverse (same
#' statement type, endpoints swapped) is dropped, so the pair is drawn once.
#' @keywords internal
#' @noRd
.consolidateEdges <- function(edges, nodes = NULL) {
    if (nrow(edges) == 0) return(edges)
    
    ptm_map      <- .ptmOverlap(edges, nodes)
    directed     <- .edge_is_directed(edges)
    consolidated <- list()
    processed    <- c()
    
    for (i in seq_len(nrow(edges))) {
        e        <- edges[i, ]
        edge_key <- paste(e$source, e$target, e$interaction, sep = "-")
        ptm_txt  <- if (edge_key %in% names(ptm_map)) ptm_map[[edge_key]] else ""
        e$category    <- .classify(e$interaction)
        e$ptm_overlap <- ptm_txt
        
        if (directed[i]) {
            e$edge_type <- "directed"
        } else {
            pair_key <- paste(c(sort(c(e$source, e$target)), e$interaction),
                              collapse = "-")
            if (pair_key %in% processed) next
            processed   <- c(processed, pair_key)
            e$edge_type <- "undirected"
        }
        consolidated[[edge_key]] <- e
    }
    
    result           <- do.call(rbind, consolidated)
    rownames(result) <- NULL
    result
}

#' Build the list of Cytoscape element objects (nodes + edges)
#'
#' Returns a list of named lists — jsonlite will serialise them cleanly.
#' Node data carry \code{status} (\code{.node_display_status()}, plus
#' \code{"sites_only"}) and \code{shape}, which the widget styles and the
#' legend lists. A protein with PTM site rows takes its colour and status from
#' its protein-level row (\code{site} empty); without one it is
#' \code{"sites_only"} and only its site nodes are coloured.
#' @keywords internal
#' @noRd
.buildElements <- function(nodes, edges, display_label_type = "id") {
    # ── node colours, status, shape ───────────────────────────────────────
    node_colors <- if ("logFC" %in% names(nodes)) {
        .mapLogFCToColor(nodes$logFC)
    } else {
        rep(NEUTRAL_NODE_COLOR, nrow(nodes))
    }
    node_status <- .node_display_status(nodes)
    node_shapes <- .node_shape(.get_column_or_default(nodes, "entity_type",
                                                     NA_character_))
    
    label_col <- if (display_label_type == "entity_name" &&
                     "entity_name" %in% names(nodes)) "entity_name" else "id"

    is_site_row <- if ("site" %in% names(nodes)) {
        !is.na(nodes$site) & trimws(nodes$site) != ""
    } else {
        rep(FALSE, nrow(nodes))
    }
    has_ptm_sites <- unique(nodes$id[is_site_row])

    # The row a protein node is drawn from: its first protein-level row,
    # else its first row
    protein_level_rows <- which(!is_site_row)
    protein_row <- protein_level_rows[match(nodes$id,
                                            nodes$id[protein_level_rows])]
    has_protein_row <- !is.na(protein_row)
    protein_row[!has_protein_row] <- match(nodes$id, nodes$id)[!has_protein_row]

    elements        <- list()
    emitted_prots   <- character(0)
    emitted_cpds    <- character(0)
    emitted_ptm_n   <- character(0)
    emitted_ptm_e   <- character(0)

    for (i in seq_len(nrow(nodes))) {
        row       <- nodes[i, , drop = FALSE]
        color     <- node_colors[i]
        has_site  <- is_site_row[i]

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
            j <- protein_row[i]
            protein_color  <- node_colors[j]
            protein_status <- node_status[j]
            if (!has_protein_row[i] && protein_status != "latent") {
                protein_color  <- NO_DATA_NODE_COLOR
                protein_status <- "sites_only"
            }
            display_label <- if (label_col == "entity_name" &&
                                 !is.na(nodes$entity_name[j]) &&
                                 nodes$entity_name[j] != "")
                nodes$entity_name[j] else row$id
            nd <- list(id        = row$id,
                       label     = display_label,
                       color     = protein_color,
                       status    = protein_status,
                       shape     = node_shapes[j],
                       entity_type = if ("entity_type" %in% names(nodes) &&
                                         !is.na(nodes$entity_type[j]))
                           nodes$entity_type[j] else "",
                       node_type = "protein",
                       width     = max(60, min(nchar(display_label) * 8 + 20, 150)),
                       height    = max(40, min(nchar(display_label) * 2 + 30, 60)))
            if (needs_compound) nd$parent <- compound_id
            elements <- c(elements, list(list(data = nd)))
            emitted_prots <- c(emitted_prots, row$id)
        }
        
        # PTM child nodes + attachment edges
        if (has_site) {
            sites <- unique(trimws(unlist(strsplit(as.character(row$site), "[_,;|]"))))
            sites <- sites[sites != ""]
            
            for (site in sites) {
                ptm_nid <- paste0(row$id, "__ptm__", site)
                if (!(ptm_nid %in% emitted_ptm_n)) {
                    elements <- c(elements, list(list(data = list(
                        id             = ptm_nid,
                        label          = site,
                        color          = color,
                        status         = node_status[i],
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
            elink <- if ("evidence_url" %in% names(row)) {
                ev <- row$evidence_url
                if (is.na(ev) || ev == "NA") "" else as.character(ev)
            } else ""
            
            elements <- c(elements, list(list(data = list(
                id          = eid,
                source      = row$source,
                target      = row$target,
                interaction = row$interaction,
                edge_type   = row$edge_type,
                category    = row$category,
                evidence_url = elink,
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

#' Resolve the displayLabelType argument of cytoscapeNetwork()
#'
#' Accepts the deprecated value "entityName" as "entity_name", with a warning.
#' @param display_label_type value passed by the caller
#' @return "id" or "entity_name"
#' @keywords internal
#' @noRd
.resolve_display_label_type <- function(display_label_type) {
    if (identical(display_label_type, "entityName")) {
        warning("displayLabelType = \"entityName\" is deprecated and will be ",
                "removed in a future release. Use \"entity_name\", which ",
                "matches the renamed nodes column.", call. = FALSE)
        return("entity_name")
    }
    if (!is.character(display_label_type) || length(display_label_type) != 1 ||
        !display_label_type %in% c("id", "entity_name")) {
        stop("displayLabelType must be \"id\" or \"entity_name\".")
    }
    display_label_type
}
