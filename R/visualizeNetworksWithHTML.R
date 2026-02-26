```r
#' Helper function to map logFC values to colors
#' @param logFC_values Numeric vector of log fold change values
#' @importFrom grDevices colorRamp rgb
#' @noRd
mapLogFCToColor <- function(logFC_values) {
    # Define the color palette
    colors <- c("#ADD8E6", "#ADD8E6", "#D3D3D3", "#FFA590", "#FFA590")
    
    # Handle case where all values are the same or missing
    if (all(is.na(logFC_values)) || length(unique(logFC_values[!is.na(logFC_values)])) <= 1) {
        return(rep("#D3D3D3", length(logFC_values)))
    }
    
    # Get range of logFC values
    default_max <- 2
    max_logFC <- max(c(abs(logFC_values), default_max), na.rm = TRUE)
    min_logFC <- -1 * max_logFC
    
    # Create color mapping function
    color_map <- colorRamp(colors)
    
    # Normalize logFC values to [0, 1] range
    normalized_values <- (logFC_values - min_logFC) / (max_logFC - min_logFC)
    
    # Handle NA values
    normalized_values[is.na(normalized_values)] <- 0.5  # Default to middle color
    
    # Get RGB colors and convert to hex
    rgb_colors <- color_map(normalized_values)
    hex_colors <- rgb(rgb_colors[,1], rgb_colors[,2], rgb_colors[,3], maxColorValue = 255)
    
    return(hex_colors)
}

# Define relationship categories and their properties
getRelationshipProperties <- function() {
    list(
        complex = list(
            types = c("Complex"),
            color = "#8B4513",        # Brown
            style = "solid",
            arrow = "none",           # Undirected
            width = 4,
            consolidate = "undirected"
        ),
        regulatory = list(
            types = c("Inhibition", "Activation", "IncreaseAmount", "DecreaseAmount"),
            colors = list(
                "Inhibition" = "#FF4444",           # Red
                "Activation" = "#44AA44",          # Green  
                "IncreaseAmount" = "#4488FF",    # Blue
                "DecreaseAmount" = "#FF8844"     # Orange
            ),
            style = "solid",
            arrow = "triangle",
            width = 3,
            consolidate = "bidirectional"
        ),
        phosphorylation = list(
            types = c("Phosphorylation"),
            color = "#9932CC",        # Purple
            style = "dashed",
            arrow = "triangle",
            width = 2,
            consolidate = "directed"
        ),
        other = list(
            color = "#666666",        # Gray
            style = "dotted",
            arrow = "triangle", 
            width = 2,
            consolidate = "directed"
        )
    )
}

#' Calculate PTM site overlap between edge targets and nodes with aggregation
#' @param edges Data frame with edge information including 'target' and 'site' columns
#' @param nodes Data frame with node information including 'id' and 'Site' columns
#' @return Vector of overlap descriptions for each unique edge (after consolidation)
#' @keywords internal
#' @noRd
.calculatePTMOverlapAggregated <- function(edges, nodes) {
    if (nrow(edges) == 0) return(character(0))
    
    # Group edges by source-target-interaction to match consolidation logic
    edges$edge_key <- paste(edges$source, edges$target, edges$interaction, sep = "-")
    unique_edges <- unique(edges$edge_key)
    
    overlap_info <- character(length(unique_edges))
    names(overlap_info) <- unique_edges
    
    for (edge_key in unique_edges) {
        # Get all edges with this source-target-interaction combination
        matching_edges <- edges[edges$edge_key == edge_key, ]
        all_overlap_sites <- c()
        
        # Process each matching edge to find PTM overlaps
        for (i in 1:nrow(matching_edges)) {
            edge <- matching_edges[i, ]
            
            # Check if edge has target and site information
            if (!is.na(edge$target) && "site" %in% names(edge) && !is.na(edge$site)) {
                # Find matching nodes with the same target ID
                target_nodes <- nodes[nodes$id == edge$target, ]
                
                if (nrow(target_nodes) > 0 && "Site" %in% names(target_nodes)) {
                    edge_sites <- trimws(unlist(strsplit(as.character(edge$site), "[,;|]")))
                    
                    # Check each target node row for site matches
                    for (j in 1:nrow(target_nodes)) {
                        if (!is.na(target_nodes$Site[j])) {
                            node_sites <- trimws(unlist(strsplit(as.character(target_nodes$Site[j]), "_")))
                            
                            # Find overlapping sites for this edge-node combination
                            overlap_sites <- intersect(edge_sites, node_sites)
                            overlap_sites <- overlap_sites[overlap_sites != "" & !is.na(overlap_sites)]
                            
                            # Add to the aggregate list
                            all_overlap_sites <- c(all_overlap_sites, overlap_sites)
                        }
                    }
                }
            }
        }
        
        # Remove duplicates and create tooltip text for this consolidated edge
        unique_overlap_sites <- unique(all_overlap_sites)
        unique_overlap_sites <- unique_overlap_sites[unique_overlap_sites != "" & !is.na(unique_overlap_sites)]
        
        if (length(unique_overlap_sites) > 0) {
            if (length(unique_overlap_sites) == 1) {
                overlap_info[edge_key] <- paste0("Overlapping PTM site: ", unique_overlap_sites[1])
            } else {
                overlap_info[edge_key] <- paste0("Overlapping PTM sites: ", paste(unique_overlap_sites, collapse = ", "))
            }
        } else {
            overlap_info[edge_key] <- ""
        }
    }
    
    return(overlap_info)
}

# Consolidate bidirectional edges based on relationship type
consolidateEdges <- function(edges, nodes = NULL) {
    if (nrow(edges) == 0) return(edges)
    
    required_cols <- c("source", "target", "interaction")
    missing_cols <- setdiff(required_cols, names(edges))
    if (length(missing_cols) > 0) {
        stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
    }
    
    # Calculate aggregated PTM overlap information if nodes are provided
    ptm_overlap_map <- if (!is.null(nodes)) {
        .calculatePTMOverlapAggregated(edges, nodes)
    } else {
        NULL
    }
    
    relationship_props <- getRelationshipProperties()
    consolidated_edges <- list()
    processed_pairs <- c()
    
    for (i in 1:nrow(edges)) {
        edge <- edges[i, ]
        pair_key <- paste(sort(c(edge$source, edge$target)), edge$interaction, collapse = "-")
        
        # Skip if we've already processed this pair
        if (pair_key %in% processed_pairs) next
        
        # Determine relationship category
        interaction_type <- edge$interaction
        category <- "other"
        for (cat_name in names(relationship_props)) {
            if (interaction_type %in% relationship_props[[cat_name]]$types) {
                category <- cat_name
                break
            }
        }
        
        # Find reverse edge if it exists
        reverse_edges <- edges[edges$source == edge$target & 
                                   edges$target == edge$source & 
                                   edges$interaction == edge$interaction, ]
        
        consolidation_type <- relationship_props[[category]]$consolidate
        
        # Get PTM overlap info for this edge combination
        edge_key <- paste(edge$source, edge$target, edge$interaction, sep = "-")
        ptm_overlap_text <- if (!is.null(ptm_overlap_map) && edge_key %in% names(ptm_overlap_map)) {
            ptm_overlap_map[[edge_key]]
        } else {
            ""
        }
        
        if (nrow(reverse_edges) > 0 && consolidation_type %in% c("undirected", "bidirectional")) {
            # Create consolidated edge
            if (consolidation_type == "undirected") {
                # For complex relationships - create undirected edge
                consolidated_edge <- data.frame(
                    source = edge$source,
                    target = edge$target,
                    interaction = edge$interaction,
                    edge_type = "undirected",
                    category = category,
                    ptm_overlap = ptm_overlap_text,
                    stringsAsFactors = FALSE
                )
            } else {
                # For regulatory relationships - create bidirectional edge
                consolidated_edge <- data.frame(
                    source = edge$source,
                    target = edge$target,
                    interaction = paste(edge$interaction, "(bidirectional)"),
                    edge_type = "bidirectional", 
                    category = category,
                    ptm_overlap = ptm_overlap_text,
                    stringsAsFactors = FALSE
                )
            }
            
            # Copy any additional columns from original edge
            other_cols <- setdiff(names(edge), c("source", "target", "interaction"))
            for (col in other_cols) {
                consolidated_edge[[col]] <- edge[[col]]
            }
            
            edge_key_final <- paste(edge$source, edge$target, consolidated_edge$interaction, sep = "-")
            consolidated_edges[[edge_key_final]] <- consolidated_edge
            
            # Mark both directions as processed
            processed_pairs <- c(processed_pairs, pair_key)
            
        } else {
            # Keep as directed edge
            directed_edge <- edge
            directed_edge$edge_type <- "directed"
            directed_edge$category <- category
            directed_edge$ptm_overlap = ptm_overlap_text
            
            edge_key_final <- paste(edge$source, edge$target, edge$interaction, sep = "-")
            consolidated_edges[[edge_key_final]] <- directed_edge
        }
    }
    
    # Convert list back to data frame
    if (length(consolidated_edges) > 0) {
        result <- do.call(rbind, consolidated_edges)
        rownames(result) <- NULL
        return(result)
    } else {
        return(edges[0, ])  # Return empty data frame with same structure
    }
}

# Get edge styling properties based on category and interaction type
getEdgeStyle <- function(interaction, category, edge_type) {
    relationship_props <- getRelationshipProperties()
    
    if (category %in% names(relationship_props)) {
        props <- relationship_props[[category]]
        
        # Handle regulatory relationships with specific colors
        if (category == "regulatory" && "colors" %in% names(props)) {
            base_interaction <- gsub(" \\(bidirectional\\)", "", interaction)
            color <- if (base_interaction %in% names(props$colors)) {
                props$colors[[base_interaction]]
            } else {
                "#666666"  # Default gray
            }
        } else {
            color <- props$color
        }
        
        # Adjust arrow type based on edge type
        arrow <- if (edge_type == "undirected") {
            "none"
        } else if (edge_type == "bidirectional") {
            "triangle"  # Will be handled specially in CSS
        } else {
            props$arrow
        }
        
        return(list(
            color = color,
            style = props$style,
            arrow = arrow,
            width = props$width
        ))
    } else {
        # Default styling for unknown relationships
        return(relationship_props$other)
    }
}

createNodeElements <- function(nodes, displayLabelType = "id") {
    if ("logFC" %in% names(nodes)) {
        node_colors <- mapLogFCToColor(nodes$logFC)
    } else {
        node_colors <- rep("#D3D3D3", nrow(nodes))
    }
    
    label_column <- if (displayLabelType == "hgncName" && "hgncName" %in% names(nodes)) {
        "hgncName"
    } else {
        "id"
    }
    
    node_elements    <- c()
    ptm_elements     <- c()
    emitted_proteins <- c()
    emitted_compounds <- c()
    emitted_ptm_nodes <- c()
    emitted_ptm_edges <- c()
    
    # Pre-compute which protein ids have at least one PTM site row,
    # so we know upfront whether a compound wrapper is needed
    has_ptm_sites <- if ("Site" %in% names(nodes)) {
        ids_with_sites <- unique(nodes$id[!is.na(nodes$Site) & trimws(nodes$Site) != ""])
        ids_with_sites
    } else {
        c()
    }
    
    for (i in seq_len(nrow(nodes))) {
        row      <- nodes[i, ]
        color    <- node_colors[i]
        has_site <- "Site" %in% names(nodes) && !is.na(row$Site) && trimws(row$Site) != ""
        
        display_label <- if (label_column == "hgncName" && !is.na(row$hgncName) && row$hgncName != "") {
            row$hgncName
        } else {
            row$id
        }
        
        needs_compound <- row$id %in% has_ptm_sites
        compound_id    <- paste0(row$id, "__compound__")
        
        # Emit invisible compound container node once per protein that has PTM children
        if (needs_compound && !(compound_id %in% emitted_compounds)) {
            node_elements <- c(node_elements,
                               paste0("{ data: { id: '", escape_js_string(compound_id),
                                      "', node_type: 'compound' } }")
            )
            emitted_compounds <- c(emitted_compounds, compound_id)
        }
        
        # Emit protein node once, assigning it to the compound if one exists
        if (!(row$id %in% emitted_proteins)) {
            parent_field <- if (needs_compound) {
                paste0(", parent: '", escape_js_string(compound_id), "'")
            } else {
                ""
            }
            node_elements <- c(node_elements,
                               paste0("{ data: { id: '", escape_js_string(row$id),
                                      "', label: '", escape_js_string(display_label),
                                      "', color: '", color,
                                      "', node_type: 'protein'",
                                      parent_field,
                                      " } }")
            )
            emitted_proteins <- c(emitted_proteins, row$id)
        }
        
        # Emit one PTM child node + attachment edge per individual site
        if (has_site) {
            sites <- trimws(unlist(strsplit(as.character(row$Site), "[_,;|]")))
            sites <- unique(sites[sites != ""])
            
            for (site in sites) {
                ptm_node_id <- paste0(row$id, "__ptm__", site)
                safe_ptm_id <- escape_js_string(ptm_node_id)
                safe_parent <- escape_js_string(row$id)
                safe_site   <- escape_js_string(site)
                
                # PTM node also belongs to the same compound container
                if (!(ptm_node_id %in% emitted_ptm_nodes)) {
                    ptm_elements <- c(ptm_elements,
                                      paste0("{ data: { id: '", safe_ptm_id,
                                             "', label: '", safe_site,
                                             "', color: '", color,
                                             "', parent_protein: '", safe_parent,
                                             "', parent: '", escape_js_string(compound_id), "'",
                                             ", node_type: 'ptm' } }")
                    )
                    emitted_ptm_nodes <- c(emitted_ptm_nodes, ptm_node_id)
                }
                
                ptm_edge_id_raw <- paste0(row$id, "__ptm_edge__", site)
                if (!(ptm_edge_id_raw %in% emitted_ptm_edges)) {
                    ptm_edge_id <- escape_js_string(ptm_edge_id_raw)
                    ptm_elements <- c(ptm_elements,
                                      paste0("{ data: { id: '", ptm_edge_id,
                                             "', source: '", safe_parent,
                                             "', target: '", safe_ptm_id,
                                             "', edge_type: 'ptm_attachment',",
                                             " category: 'ptm_attachment',",
                                             " interaction: '',",
                                             " color: '", color, "',",
                                             " line_style: 'dotted',",
                                             " arrow_shape: 'none',",
                                             " width: 1.5,",
                                             " tooltip: '' } }")
                    )
                    emitted_ptm_edges <- c(emitted_ptm_edges, ptm_edge_id_raw)
                }
            }
        }
    }
    
    return(c(node_elements, ptm_elements))
}

createEdgeElements <- function(edges, nodes = NULL) {
    if (nrow(edges) == 0) return(list())
    
    # First consolidate edges
    consolidated_edges <- consolidateEdges(edges, nodes)
    
    edge_elements <- list()
    
    for (i in 1:nrow(consolidated_edges)) {
        row <- consolidated_edges[i,]
        edge_key <- paste(row$source, row$target, row$interaction, sep = "-")
        
        # Get styling for this edge
        style <- getEdgeStyle(row$interaction, row$category, row$edge_type)
        
        # Sanitize optional evidenceLink
        evidence_link <- if ("evidenceLink" %in% names(row)) row$evidenceLink else NA_character_
        evidence_link <- ifelse(is.na(evidence_link) | evidence_link == "NA", "", evidence_link)
        evidence_link <- escape_js_string(evidence_link)
        
        # Escape quotes in tooltip text for JavaScript safety
        tooltip_text <- gsub("'", "\\\\'", row$ptm_overlap)
        
        # Create edge data with styling information and PTM overlap tooltip
        edge_data <- paste0("{ data: { source: '", row$source, 
                            "', target: '", row$target, 
                            "', id: '", edge_key,
                            "', interaction: '", row$interaction,
                            "', edge_type: '", row$edge_type,
                            "', category: '", row$category,
                            "', evidenceLink: '", evidence_link,
                            "', color: '", style$color,
                            "', line_style: '", style$style,
                            "', arrow_shape: '", style$arrow,
                            "', width: ", style$width,
                            ", tooltip: '", tooltip_text, "' } }")
        
        edge_elements[[edge_key]] <- edge_data
    }
    
    return(edge_elements)
}

# Helper to safely embed strings into JS single-quoted literals
escape_js_string <- function(x) {
    if (is.null(x)) return("")
    x <- as.character(x)
    x <- gsub("\\\\", "\\\\\\\\", x)  # backslashes
    x <- gsub("'", "\\\\'", x)        # single quotes
    x <- gsub("\r", "\\\\r", x)
    x <- gsub("\n", "\\\\n", x)
    x
}

#' Generate Cytoscape visualization configuration
#' 
#' This function creates a complete Cytoscape configuration object that can be
#' used to render a network visualization. It's decoupled from any specific
#' UI framework.
#' 
#' @param nodes List of nodes from getSubnetworkFromIndra
#' @param edges List of edges from getSubnetworkFromIndra
#' @param display_label_type column of nodes table for displaying node names
#' @param container_id ID of the HTML container element (default: 'network-cy')
#' @param event_handlers Optional list of event handler configurations
#' @param layout_options Optional list of layout configuration options
#' @export
#' @return List containing:
#'   - elements: Combined node and edge elements
#'   - style: Cytoscape style configuration
#'   - layout: Layout configuration
#'   - container_id: Container element ID
#'   - js_code: Complete JavaScript code (for backward compatibility)
generateCytoscapeConfig <- function(nodes, edges, 
                                    display_label_type = "id",
                                    container_id = "network-cy",
                                    event_handlers = NULL,
                                    layout_options = NULL,
                                    node_font_size = 12) {
    
    # Create elements
    node_elements <- createNodeElements(nodes, display_label_type)
    edge_elements <- createEdgeElements(edges, nodes)
    
    # Default layout options
    default_layout <- list(
        name = "dagre",
        rankDir = "TB",
        animate = TRUE,
        fit = TRUE,
        padding = 30,
        spacingFactor = 1.5,
        nodeSep = 50,
        edgeSep = 20,
        rankSep = 80
    )
    
    # Merge with custom layout options if provided
    layout_config <- default_layout
    if (!is.null(layout_options)) {
        for (layout_name in names(layout_options)) {
            layout_config[[layout_name]] <- layout_options[[layout_name]]
        }
    }
    
    # Define the style configuration
    style_config <- list(
        list(
            selector = "node",
            style = list(
                `background-color` = "data(color)",
                label = "data(label)",
                width = "function(ele) { var label = ele.data('label') || ''; var labelLength = label.length; return Math.max(60, Math.min(labelLength * 8 + 20, 150)); }",
                height = "function(ele) { var label = ele.data('label') || ''; var labelLength = label.length; return Math.max(40, Math.min(labelLength * 2 + 30, 60)); }",
                shape = "round-rectangle",
                `font-size` = paste0(node_font_size, "px"),
                `font-weight` = "bold",
                color = "#000",
                `text-valign` = "center",
                `text-halign` = "center",
                `text-wrap` = "wrap",
                `text-max-width` = "function(ele) { var label = ele.data('label') || ''; var labelLength = label.length; return Math.max(50, Math.min(labelLength * 8 + 10, 140)); }",
                `border-width` = 2,
                `border-color` = "#333",
                padding = "5px"
            )
        ),
        list(
            selector = "edge",
            style = list(
                width = "data(width)",
                `line-color` = "data(color)",
                `line-style` = "data(line_style)",
                label = "data(interaction)",
                `curve-style` = "bezier",
                `target-arrow-shape` = "data(arrow_shape)",
                `target-arrow-color` = "data(color)",
                `source-arrow-shape` = "function(ele) { return ele.data('edge_type') === 'bidirectional' ? 'triangle' : 'none'; }",
                `source-arrow-color` = "data(color)",
                `edge-text-rotation` = "autorotate",
                `text-margin-y` = -12,
                `text-halign` = "center",
                `font-size` = "16px",
                `font-weight` = "bold",
                color = "data(color)",
                `text-background-color` = "#ffffff",
                `text-background-opacity` = 0.8,
                `text-background-padding` = "2px"
            )
        ),
        list(
            selector = "edge[category = 'complex']",
            style = list(
                `line-style` = "solid",
                `target-arrow-shape` = "none",
                `source-arrow-shape` = "none"
            )
        ),
        list(
            selector = "edge[category = 'phosphorylation']",
            style = list(
                `line-style` = "dashed",
                width = 2
            )
        ),
        list(
            selector = "edge[edge_type = 'bidirectional']",
            style = list(
                `source-arrow-shape` = "triangle",
                `target-arrow-shape` = "triangle"
            )
        ),
        list(
            selector = "node[node_type = 'ptm']",
            style = list(
                shape = "ellipse",
                width = "20px",
                height = "20px",
                `background-color` = "data(color)",   # <-- was hardcoded "#9932CC"
                `border-color` = "#333",              # <-- neutral border instead of purple
                `border-width` = 1.5,
                label = "data(label)",
                `font-size` = "8px",
                `font-weight` = "normal",
                color = "#000000",                    # <-- dark text works across the logFC palette
                `text-valign` = "center",
                `text-halign` = "center",
                `text-wrap` = "wrap",
                `text-max-width` = "18px"
            )
        ),
        # --- PTM attachment edge style (hide label, keep it subtle) ---
        list(
            selector = "edge[edge_type = 'ptm_attachment']",
            style = list(
                `line-style` = "dotted",
                `line-color` = "#9932CC",
                width = 1.5,
                `target-arrow-shape` = "none",
                `source-arrow-shape` = "none",
                label = "",             # no label on these connector edges
                `z-index` = 0          # render behind main edges
            )
        ),
        list(
            selector = "node[node_type = 'compound']",
            style = list(
                `background-opacity` = 0,
                `border-width` = 0,
                `border-opacity` = 0,
                `padding` = "10px",
                label = "",
                `z-index` = 0
            )
        )
    )
    
    # Combine elements
    elements <- c(node_elements, edge_elements)
    
    # Create the main configuration object
    config <- list(
        elements = elements,
        style = style_config,
        layout = layout_config,
        container_id = container_id,
        event_handlers = event_handlers
    )
    
    # Generate JavaScript code for backward compatibility
    config$js_code <- generateJavaScriptCode(config)
    
    return(config)
}

#' Generate JavaScript code from Cytoscape configuration
#' 
#' Internal function to convert configuration object to JavaScript code
#' 
#' @param config Configuration object from generateCytoscapeConfig()
#' @return Character string containing JavaScript code
generateJavaScriptCode <- function(config) {
    
    # Convert R list to JSON-like string for JavaScript
    elements_js <- paste(config$elements, collapse = ", ")
    
    # Convert style configuration to JavaScript
    style_js <- convertStyleToJS(config$style)
    
    # Convert layout configuration to JavaScript
    layout_js <- convertLayoutToJS(config$layout)
    
    # Build event handlers JavaScript
    event_handlers_js <- ""
    if (!is.null(config$event_handlers)) {
        handlers <- sapply(names(config$event_handlers), function(event) {
            handler_code <- config$event_handlers[[event]]
            switch(event,
                   "edge_click" = paste0("cy.on('tap', 'edge', ", handler_code, ");"),
                   "node_click" = paste0("cy.on('tap', 'node', ", handler_code, ");"),
                   handler_code  # Custom event handler
            )
        })
        event_handlers_js <- paste(handlers, collapse = "\n    ")
    }
    
    # Generate the complete JavaScript code
    js_code <- paste0("
    cytoscape.use(cytoscapeDagre);
    var cy = cytoscape({
        container: document.getElementById('", config$container_id, "'),
        elements: [", elements_js, "],
        style: ", style_js, ",
        layout: ", layout_js, "
    });
    
    // After layout completes, reposition PTM nodes directly beside their parent protein
    cy.on('layoutstop', function() {
        var ptmNodes = cy.nodes('[node_type = \"ptm\"]');
        ptmNodes.forEach(function(ptmNode) {
            var parentId = ptmNode.data('parent_protein');
            var parentNode = cy.getElementById(parentId);
            if (parentNode.length === 0) return;

            var parentPos = parentNode.position();
            var parentW   = parentNode.outerWidth();
            var parentH   = parentNode.outerHeight();
            var ptmR      = ptmNode.outerWidth() / 2;  // PTM node is a small circle

            // Collect all PTM siblings so we can fan them around the parent
            var siblings = cy.nodes('[parent_protein = \"' + parentId + '\"]');
            var idx      = siblings.indexOf(ptmNode);
            var total    = siblings.length;

            // Distribute siblings evenly across the bottom arc of the parent
            // angleStart/End in radians: spread across bottom 180 degrees
            var angleStart = Math.PI * 0.15;
            var angleEnd   = Math.PI * 0.85;
            var angle = total === 1
                ? Math.PI / 2   // single PTM: directly below center
                : angleStart + (angleEnd - angleStart) * (idx / (total - 1));

            // Place PTM node just outside the parent border
            var offsetX = (parentW / 2 + ptmR + 4) * Math.cos(angle);
            var offsetY = (parentH / 2 + ptmR + 4) * Math.sin(angle);

            ptmNode.position({
                x: parentPos.x + offsetX,
                y: parentPos.y + offsetY
            });
        });
    });
    
    // Create tooltip element
    var tooltip = document.createElement('div');
    tooltip.style.cssText = `
        position: absolute;
        background-color: rgba(0, 0, 0, 0.9);
        color: white;
        padding: 8px 12px;
        border-radius: 4px;
        font-size: 12px;
        font-family: Arial, sans-serif;
        white-space: nowrap;
        pointer-events: none;
        z-index: 9999;
        box-shadow: 0 2px 8px rgba(0, 0, 0, 0.3);
        display: none;
        max-width: 300px;
        word-wrap: break-word;
        white-space: pre-wrap;
    `;
    document.body.appendChild(tooltip);
    
    // Only show tooltip if there's actual PTM overlap information
    cy.on('mouseover', 'edge', function(evt) {
        var edge = evt.target;
        var tooltipText = edge.data('tooltip');
        if (tooltipText && tooltipText.trim() !== '' && tooltipText.trim() !== 'No overlapping PTM sites found') {
            tooltip.innerHTML = tooltipText;
            tooltip.style.display = 'block';
        }
    });
    
    cy.on('mousemove', 'edge', function(evt) {
        if (tooltip.style.display === 'block') {
            tooltip.style.left = evt.originalEvent.pageX + 10 + 'px';
            tooltip.style.top = evt.originalEvent.pageY - 30 + 'px';
        }
    });
    
    cy.on('mouseout', 'edge', function(evt) {
        tooltip.style.display = 'none';
    });
    
    ", event_handlers_js)
    
    return(js_code)
}

# Helper function to convert style list to JavaScript
convertStyleToJS <- function(style_list) {
    style_items <- sapply(style_list, function(item) {
        # Properly escape selector strings, especially those with special characters
        selector_js <- paste0("\"", gsub("\"", "\\\"", item$selector), "\"")
        
        # Convert style properties
        style_props <- sapply(names(item$style), function(prop) {
            value <- item$style[[prop]]
            if (is.character(value) && !grepl("^function\\(", value)) {
                # Use double quotes and escape any existing double quotes
                escaped_prop <- gsub("\"", "\\\"", prop)
                escaped_value <- gsub("\"", "\\\"", value)
                paste0("\"", escaped_prop, "\": \"", escaped_value, "\"")
            } else {
                escaped_prop <- gsub("\"", "\\\"", prop)
                paste0("\"", escaped_prop, "\": ", value)
            }
        })
        
        paste0("{ selector: ", selector_js, ", style: { ", paste(style_props, collapse = ", "), " } }")
    })
    
    paste0("[", paste(style_items, collapse = ", "), "]")
}

# Helper function to convert layout list to JavaScript
convertLayoutToJS <- function(layout_list) {
    layout_props <- sapply(names(layout_list), function(prop) {
        value <- layout_list[[prop]]
        if (is.character(value)) {
            escaped_prop <- gsub("\"", "\\\"", prop)
            escaped_value <- gsub("\"", "\\\"", value)
            paste0("\"", escaped_prop, "\": \"", escaped_value, "\"")
        } else if (is.logical(value)) {
            escaped_prop <- gsub("\"", "\\\"", prop)
            paste0("\"", escaped_prop, "\": ", tolower(value))
        } else {
            escaped_prop <- gsub("\"", "\\\"", prop)
            paste0("\"", escaped_prop, "\": ", value)
        }
    })
    
    paste0("{ ", paste(layout_props, collapse = ", "), " }")
}

#' Export network data with Cytoscape visualization
#' 
#' Convenience function that takes nodes and edges data directly and creates
#' both the configuration and HTML export in one step.
#' 
#' @param nodes \code{data.frame} with node information.
#' @param edges \code{data.frame} with edge information.
#' @param filename \code{character} specifying the output HTML filename.
#' @param displayLabelType \code{character} specifying the type of label to display ("id" or "hgncName").
#' @param nodeFontSize \code{numeric} specifying the font size for node labels (default: 12).
#' @param ... Additional arguments passed to \code{exportCytoscapeToHTML()}.
#' @return Invisibly returns the file path of the created HTML file.
#' @export
#' @examples
#' \dontrun{
#' # Example nodes and edges data frames
#' nodes <- data.frame(id = c("A", "B", "C"), logFC = c(1.5, -0.5, 0))
#' edges <- data.frame(source = c("A", "B"), target = c("B", "C"), interaction = c("Activation", "Inhibition"))
#' 
#' # Export network to HTML
#' exportNetworkToHTML(nodes, edges, filename = "my_network.html")
#' }
exportNetworkToHTML <- function(nodes, edges, 
                                filename = "network_visualization.html",
                                displayLabelType = "id",
                                nodeFontSize = 12,
                                ...) {
    
    event_handlers <- createEdgeClickHandler()
    
    # Generate configuration
    config <- generateCytoscapeConfig(
        nodes, 
        edges, 
        display_label_type = displayLabelType, 
        node_font_size = nodeFontSize,
        event_handlers = event_handlers
    )
    
    # Export to HTML
    exportCytoscapeToHTML(config, filename, ...)
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
    
    event_handlers <- createEdgeClickHandler()
    
    # Generate configuration
    config <- generateCytoscapeConfig(
        nodes, 
        edges, 
        display_label_type = displayLabelType, 
        node_font_size = nodeFontSize,
        event_handlers = event_handlers
    )
    
    # Create temporary filename
    temp_file <- tempfile(fileext = ".html")
    
    # Export to temp file
    exportCytoscapeToHTML(config, temp_file, ...)
    
    # Open in browser
    if (interactive()) {
        browseURL(temp_file)
        cat("Network opened in browser. Temporary file:", temp_file, "\n")
    }
    
    invisible(temp_file)
}
```