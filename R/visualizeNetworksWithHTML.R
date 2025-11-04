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
#' @noRd
calculatePTMOverlapAggregated <- function(edges, nodes) {
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
        calculatePTMOverlapAggregated(edges, nodes)
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
    # Map logFC to colors if logFC column exists
    if ("logFC" %in% names(nodes)) {
        node_colors <- mapLogFCToColor(nodes$logFC)
    } else {
        node_colors <- rep("#D3D3D3", nrow(nodes))  # Default color
    }
    
    # Determine which column to use for labels
    label_column <- if(displayLabelType == "hgncName" && "hgncName" %in% names(nodes)) {
        "hgncName"
    } else {
        "id"
    }
    
    apply(cbind(nodes, color = node_colors), 1, function(row) {
        # Use the appropriate label, fallback to id if hgncName is missing/empty
        display_label <- if(label_column == "hgncName" && !is.na(row['hgncName']) && row['hgncName'] != "") {
            row['hgncName']
        } else {
            row['id']
        }
        
        paste0("{ data: { id: '", row['id'], "', label: '", display_label, "', color: '", row['color'], "' } }")
    })
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

#' Export Cytoscape network visualization to standalone HTML file
#' 
#' This function takes a Cytoscape configuration object and creates a complete
#' standalone HTML file that can be opened in any web browser.
#' 
#' @param config Configuration object from generateCytoscapeConfig()
#' @param filename Output HTML filename (default: "network_visualization.html")
#' @param title HTML page title (default: "Network Visualization")
#' @param width Container width (default: "100%")
#' @param height Container height (default: "600px")
#' @param include_controls Whether to include basic zoom/fit controls (default: TRUE)
#' @param custom_css Additional CSS styling (optional)
#' @param custom_js Additional JavaScript code (optional)
#' 
#' @return Invisibly returns the file path of the created HTML file
#' 
#' @examples
#' \dontrun{
#' # Assuming you have nodes and edges data
#' config <- generateCytoscapeConfig(node_elements, edge_elements)
#' 
#' # Export to HTML
#' exportCytoscapeToHTML(config, "my_network.html")
#' }
#' @noRd
exportCytoscapeToHTML <- function(config, 
                                  filename = "network_visualization.html",
                                  title = "Network Visualization",
                                  width = "100%", 
                                  height = "600px",
                                  include_controls = TRUE,
                                  custom_css = "",
                                  custom_js = "") {
    
    # Validate config object
    if (!is.list(config) || !all(c("elements", "style", "layout", "container_id") %in% names(config))) {
        stop("Invalid config object. Must be generated by generateCytoscapeConfig()")
    }
    
    # Generate the JavaScript code if not already present
    if (!"js_code" %in% names(config)) {
        config$js_code <- generateJavaScriptCode(config)
    }
    
    # Create controls HTML and JavaScript if requested
    controls_html <- ""
    controls_js <- ""
    controls_css <- ""
    
    if (include_controls) {
        controls_html <- '
    <div id="controls" style="margin-bottom: 10px;">
      <button id="fit-btn" class="control-btn">Fit to Screen</button>
      <button id="center-btn" class="control-btn">Center</button>
      <button id="zoom-in-btn" class="control-btn">Zoom In</button>
      <button id="zoom-out-btn" class="control-btn">Zoom Out</button>
      <button id="reset-btn" class="control-btn">Reset View</button>
      <button id="export-png-btn" class="control-btn export-btn">Export as PNG</button>
    </div>'
        
        controls_css <- '
    .control-btn {
      margin: 2px;
      padding: 6px 12px;
      background-color: #f8f9fa;
      border: 1px solid #dee2e6;
      border-radius: 4px;
      cursor: pointer;
      font-size: 12px;
    }
    .control-btn:hover {
      background-color: #e9ecef;
    }
    .control-btn:active {
      background-color: #dee2e6;
    }
    .export-btn {
      background-color: #28a745;
      color: white;
      border-color: #28a745;
      font-weight: bold;
    }
    .export-btn:hover {
      background-color: #218838;
      border-color: #1e7e34;
    }
    .export-btn:active {
      background-color: #1e7e34;
    }'
        
        controls_js <- '
    // Add control event listeners
    document.getElementById("fit-btn").addEventListener("click", function() {
      cy.fit();
    });
    
    document.getElementById("center-btn").addEventListener("click", function() {
      cy.center();
    });
    
    document.getElementById("zoom-in-btn").addEventListener("click", function() {
      cy.zoom({
        level: cy.zoom() * 1.2,
        renderedPosition: { x: cy.width()/2, y: cy.height()/2 }
      });
    });
    
    document.getElementById("zoom-out-btn").addEventListener("click", function() {
      cy.zoom({
        level: cy.zoom() * 0.8,
        renderedPosition: { x: cy.width()/2, y: cy.height()/2 }
      });
    });
    
    document.getElementById("reset-btn").addEventListener("click", function() {
      cy.reset();
      cy.fit();
    });
    
    // Export PNG functionality
    document.getElementById("export-png-btn").addEventListener("click", function() {
      const exportBtn = this;
      const originalText = exportBtn.textContent;
      
      // Disable button and show loading state
      exportBtn.disabled = true;
      exportBtn.textContent = "Exporting...";
      
      try {
        // Generate high-resolution PNG (scale factor of 8 for high quality)
        const png64 = cy.png({
          output: "base64uri",
          bg: "white",
          full: true,  // Export the entire graph
          scale: 10,    // 10x resolution for publication quality
          maxWidth: 10000,  // Maximum width in pixels
          maxHeight: 10000  // Maximum height in pixels
        });
        
        // Create download link
        const link = document.createElement("a");
        link.href = png64;
        link.download = "network_visualization_" + new Date().getTime() + ".png";
        
        // Trigger download
        document.body.appendChild(link);
        link.click();
        document.body.removeChild(link);
        
        // Show success feedback
        exportBtn.textContent = "✓ Exported!";
        exportBtn.style.backgroundColor = "#28a745";
        
        // Reset button after 2 seconds
        setTimeout(function() {
          exportBtn.disabled = false;
          exportBtn.textContent = originalText;
          exportBtn.style.backgroundColor = "";
        }, 2000);
        
      } catch (error) {
        console.error("Error exporting PNG:", error);
        alert("Error exporting PNG: " + error.message);
        exportBtn.disabled = false;
        exportBtn.textContent = originalText;
      }
    });'
    }
    
    # Create the complete HTML content
    html_content <- paste0('<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>', title, '</title>
    
    <!-- Cytoscape.js and dependencies -->
    <script src="https://cdnjs.cloudflare.com/ajax/libs/cytoscape/3.32.0/cytoscape.min.js"></script>
    <script src="https://cdnjs.cloudflare.com/ajax/libs/graphlib/2.1.8/graphlib.min.js"></script>
    <script src="https://cdnjs.cloudflare.com/ajax/libs/dagre/0.8.5/dagre.min.js"></script>
    <script src="https://unpkg.com/cytoscape-dagre@2.3.0/cytoscape-dagre.js"></script>
    
    <style>
        body {
            font-family: Arial, sans-serif;
            margin: 0;
            padding: 20px;
            background-color: #f5f5f5;
        }
        
        .container {
            max-width: 1200px;
            margin: 0 auto;
            background-color: white;
            border-radius: 8px;
            box-shadow: 0 2px 10px rgba(0,0,0,0.1);
            padding: 20px;
        }
        
        h1 {
            color: #333;
            text-align: center;
            margin-bottom: 20px;
        }
        
        #'
                           , config$container_id, ' {
            width: ', width, ';
            height: ', height, ';
            border: 1px solid #ddd;
            border-radius: 4px;
            background-color: #fff;
        }
        
        #legend {
            background-color: #f8f9fa;
            border: 1px solid #dee2e6;
            border-radius: 4px;
            padding: 15px;
        }
        
        .legend-title {
            font-weight: bold;
            margin-bottom: 10px;
            font-size: 16px;
            color: #333;
        }
        
        .legend-item {
            display: flex;
            align-items: center;
            margin-bottom: 8px;
            font-size: 12px;
        }
        
        .legend-color {
            width: 20px;
            height: 20px;
            border: 2px solid #333;
            border-radius: 3px;
            margin-right: 8px;
        }
        
        .legend-gradient {
            height: 120px;
            width: 20px;
            border: 2px solid #333;
            border-radius: 3px;
            margin-right: 8px;
            background: linear-gradient(to top, #ADD8E6, #D3D3D3, #FFA590);
        }
        
        .legend-gradient-labels {
            display: flex;
            flex-direction: column;
            justify-content: space-between;
            height: 120px;
            font-size: 14px;
        }
        
        .edge-legend {
            margin-top: 20px;
        }
        
        .edge-legend-item {
            display: flex;
            align-items: center;
            margin-bottom: 6px;
            font-size: 14px;
        }
        
        .edge-legend-line {
            width: 30px;
            height: 2px;
            margin-right: 8px;
        }
        
        ', controls_css, '
        
        ', custom_css, '
        
        .info-panel {
            margin-top: 15px;
            padding: 10px;
            background-color: #f8f9fa;
            border: 1px solid #dee2e6;
            border-radius: 4px;
            font-size: 12px;
            color: #666;
        }
    </style>
</head>
<body>
    <div class="container">
        <h1>', title, '</h1>
        
        ', controls_html, '
        
        <div style="display: flex; gap: 20px;">
            <div id="', config$container_id, '" style="flex: 1;"></div>
            <div id="legend" style="width: 200px; flex-shrink: 0;">
                <!-- Legend will be populated by JavaScript -->
            </div>
        </div>
        
        <div class="info-panel">
            <strong>Instructions:</strong> 
            Click and drag to pan the network
            | Use mouse wheel to zoom in/out
            | <strong>Hover over edges to see PTM site overlap information</strong>
            | Click on nodes or edges to select them
            ', if(include_controls) '| Use the buttons above for common navigation actions' else '', '
            ', if(include_controls) '| <strong>Export as PNG to save the current view in high resolution</strong>' else '', '
        </div>
    </div>
    
    <script>
        // Function to create the legend
        function createLegend(cy) {
            const legendDiv = document.getElementById("legend");
            
            // Check if nodes have logFC data
            const nodes = cy.nodes();
            let hasLogFC = false;
            let logFCValues = [];
            
            nodes.forEach(function(node) {
                const nodeData = node.data();
                // Try to extract logFC from color or check if we can determine logFC values
                // Since we only have access to the final colors, we\'ll create a standard legend
                hasLogFC = true; // Assume we have logFC if we\'re showing the legend
            });
            
            let legendHTML = "";
            
            if (hasLogFC) {
                legendHTML += `
                    <div class="legend-title">Node Colors (logFC)</div>
                    <div class="legend-item">
                        <div class="legend-gradient"></div>
                        <div class="legend-gradient-labels">
                            <div>Upregulated</div>
                            <div>Neutral (0)</div>
                            <div>Downregulated</div>
                        </div>
                    </div>
                    <div style="margin-top: 10px; font-size: 14px; color: #666;">
                        Log Fold Change values
                    </div>
                `;
            } else {
                legendHTML += `
                    <div class="legend-title">Node Colors</div>
                    <div class="legend-item">
                        <div class="legend-color" style="background-color: #D3D3D3;"></div>
                        <span>Default</span>
                    </div>
                `;
            }
            
            const edges = cy.edges();
            
            const edgeTypeConfigs = [
                { type: "Activation", color: "#44AA44", label: "Activation", style: "" },
                { type: "Inhibition", color: "#FF4444", label: "Inhibition", style: "" },
                { type: "IncreaseAmount", color: "#4488FF", label: "Increase Amount", style: "" },
                { type: "DecreaseAmount", color: "#FF8844", label: "Decrease Amount", style: "" },
                { type: "Phosphorylation", color: "#9932CC", label: "Phosphorylation", style: "border-style: dashed; border-width: 1px; height: 0px; border-top-width: 2px;" },
                { type: "Complex", color: "#8B4513", label: "Complex", style: "" }
            ];
            
            const existingEdgeTypes = new Set();
            edges.forEach(edge => {
                const edgeType = edge.data("interaction");
                if (edgeType) {
                    const normalizedType = edgeType.replace(" (bidirectional)", "");
                    existingEdgeTypes.add(normalizedType);
                }
            });
            
            let edgeLegendItems = "";
            edgeTypeConfigs.forEach(config => {
                if (existingEdgeTypes.has(config.type)) {
                    const styleAttr = config.style ? `style="background-color: ${config.color}; ${config.style}"` : `style="background-color: ${config.color};"`;
                    edgeLegendItems += `
                                <div class="edge-legend-item">
                                    <div class="edge-legend-line" ${styleAttr}></div>
                                    <span>${config.label}</span>
                                </div>`;
                }
            });
            
            if (edgeLegendItems) {
                legendHTML += `
                    <div class="edge-legend">
                        <div class="legend-title">Edge Types</div>${edgeLegendItems}
                    </div>
                    <div style="margin-top: 15px; padding: 8px; background-color: #e3f2fd; border-radius: 4px; font-size: 10px;">
                        <strong>PTM Site Info:</strong><br>
                        Hover over edges to see overlapping PTM sites between the edge target and node data
                    </div>
                `;
            }
            
            legendDiv.innerHTML = legendHTML;
        }
        
        // Wait for DOM to be fully loaded
        document.addEventListener("DOMContentLoaded", function() {
            try {
                // Initialize Cytoscape
                ', config$js_code, '
                
                // Add basic interactivity
                cy.on("tap", "node", function(evt) {
                    var node = evt.target;
                    console.log("Node clicked:", node.data());
                });
                
                cy.on("tap", "edge", function(evt) {
                    var edge = evt.target;
                    console.log("Edge clicked:", edge.data());
                });
                
                ', controls_js, '
                
                // Create legend
                createLegend(cy);
                
                ', custom_js, '
                
                console.log("Network visualization loaded successfully!");
                
            } catch (error) {
                console.error("Error loading network visualization:", error);
                document.getElementById("', config$container_id, '").innerHTML = 
                    "<div style=\\"padding: 20px; text-align: center; color: red;\\">"+
                    "Error loading visualization: " + error.message + "</div>";
            }
        });
    </script>
</body>
</html>')
    
    # Write the HTML file
    writeLines(html_content, filename)
    
    # Print success message
    cat("Network visualization exported to:", normalizePath(filename), "\n")
    
    # Return the file path invisibly
    invisible(normalizePath(filename))
}

#' Export network data with Cytoscape visualization
#' 
#' Convenience function that takes nodes and edges data directly and creates
#' both the configuration and HTML export in one step.
#' 
#' @param nodes Data frame with node information
#' @param edges Data frame with edge information  
#' @param filename Output HTML filename
#' @param displayLabelType Type of label to display ("id" or "hgncName")
#' @param ... Additional arguments passed to exportCytoscapeToHTML()
#' @export
#' @return Invisibly returns the file path of the created HTML file
exportNetworkToHTML <- function(nodes, edges, 
                                filename = "network_visualization.html",
                                displayLabelType = "id",
                                nodeFontSize = 12,
                                ...) {
    
    # Generate configuration
    config <- generateCytoscapeConfig(nodes, edges, display_label_type = displayLabelType, node_font_size = nodeFontSize)
    
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
    
    # Generate configuration
    config <- generateCytoscapeConfig(nodes, edges, display_label_type = displayLabelType, node_font_size = nodeFontSize)
    
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