#' Map logFC values to a blue-grey-red colour palette
#' @importFrom grDevices colorRamp rgb
#' @keywords internal
#' @noRd
.mapLogFCToColor <- function(logFC_values) {
    colors <- c("#ADD8E6", "#ADD8E6", "#D3D3D3", "#FFA590", "#FFA590")
    
    if (all(is.na(logFC_values)) ||
        length(unique(logFC_values[!is.na(logFC_values)])) <= 1) {
        return(rep("#D3D3D3", length(logFC_values)))
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
            consolidate = "directed"
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
    
    color <- if (category == "regulatory" && !is.null(p$colors)) {
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

#' Aggregate PTM overlap between edge targets and node Site columns
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

#' Consolidate undirected edges
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
        
        if (nrow(rev_edges) > 0 && con_type == "undirected") {
            new_edge <- data.frame(source      = e$source,
                                   target      = e$target,
                                   interaction = e$interaction,
                                   edge_type   = "undirected",
                                   category    = cat,
                                   ptm_overlap = ptm_txt,
                                   stringsAsFactors = FALSE)
            for (col in setdiff(names(e), c("source", "target", "interaction"))) {
                new_edge[[col]] <- e[[col]]
            }
            key <- paste(e$source, e$target, e$interaction, sep = "-")
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
    
    label_col <- if (display_label_type == "entityName" &&
                     "entityName" %in% names(nodes)) "entityName" else "id"

    has_ptm_sites <- if ("Site" %in% names(nodes)) {
        unique(nodes$id[!is.na(nodes$Site) & trimws(nodes$Site) != ""])
    } else {
        character(0)
    }

    elements        <- list()
    emitted_prots   <- character(0)
    # `emitted_cpds` and `node_type = "compound"` below refer to Cytoscape
    # grouping containers used to parent PTM satellite nodes around a protein.
    # This Cytoscape "compound" concept is UNRELATED to the chemical
    # `proteinIdType = "Compound"` analyte type in annotateProteinInfoFromIndra.
    emitted_cpds    <- character(0)
    emitted_ptm_n   <- character(0)
    emitted_ptm_e   <- character(0)

    for (i in seq_len(nrow(nodes))) {
        row       <- nodes[i, , drop = FALSE]
        color     <- node_colors[i]
        has_site  <- "Site" %in% names(nodes) &&
            !is.na(row$Site) && trimws(row$Site) != ""

        display_label <- if (label_col == "entityName" &&
                             !is.na(row$entityName) && row$entityName != "")
            row$entityName else row$id

        needs_compound <- row$id %in% has_ptm_sites
        compound_id    <- paste0(row$id, "__compound__")

        # Cytoscape compound container (PTM grouping parent — not a chemical compound)
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
                       node_type = "protein",
                       width     = max(60, min(nchar(display_label) * 8 + 20, 150)),
                       height    = max(40, min(nchar(display_label) * 2 + 30, 60)))
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