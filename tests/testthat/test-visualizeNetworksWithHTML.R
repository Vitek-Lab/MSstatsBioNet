# =============================================================================
# UNIT TESTS FOR NETWORK VISUALIZATION MODULE
# =============================================================================

# Load required libraries

# =============================================================================
# MOCK DATA SETUP
# =============================================================================

# Mock data for testing
create_mock_input_data <- function() {
    data.frame(
        Protein = c("P53_HUMAN", "MDM2_HUMAN", "ATM_HUMAN", "BRCA1_HUMAN"),
        log2FC = c(2.5, -1.8, 1.2, -2.1),
        adj.pvalue = c(0.001, 0.02, 0.03, 0.005),
        Label = rep("Treatment_vs_Control", 4),
        stringsAsFactors = FALSE
    )
}

create_mock_annotated_data <- function() {
    data.frame(
        Protein = c("P53_HUMAN", "MDM2_HUMAN", "ATM_HUMAN", "BRCA1_HUMAN"),
        log2FC = c(2.5, -1.8, 1.2, -2.1),
        adj.pvalue = c(0.001, 0.02, 0.03, 0.005),
        Label = rep("Treatment_vs_Control", 4),
        HgncId = c("101", "102", "103", "104"),
        HgncName = c("TP53", "MDM2", "ATM", "BRCA1"),
        stringsAsFactors = FALSE
    )
}

create_mock_subnetwork_nodes <- function() {
    data.frame(
        id = c("P53_HUMAN", "MDM2_HUMAN", "ATM_HUMAN", "BRCA1_HUMAN"),
        logFC = c(2.5, -1.8, 1.2, -2.1),
        pvalue = c(0.001, 0.02, 0.03, 0.005),
        hgncName = c("TP53", "MDM2", "ATM", "BRCA1"),
        stringsAsFactors = FALSE
    )
}

create_mock_subnetwork_edges <- function() {
    data.frame(
        source = c("TP53", "MDM2", "ATM", "TP53"),
        target = c("MDM2", "TP53", "TP53", "BRCA1"),
        interaction = c("Inhibition", "Activation", "Phosphorylation", "Complex"),
        evidenceCount = c(15, 8, 12, 5),
        evidenceLink = c("link1", "link2", "link3", "link4"),
        source_counts = c("{reach:10, signor:5}", "{reach:5,biopax:3}", "{reach:8,phosphoelm:4}", "{biopax:5}"),
        stringsAsFactors = FALSE
    )
}

create_mock_subnetwork <- function() {
    list(
        nodes = create_mock_subnetwork_nodes(),
        edges = create_mock_subnetwork_edges()
    )
}

# =============================================================================
# TESTS FOR COLOR MAPPING FUNCTION
# =============================================================================

test_that("mapLogFCToColor handles various input scenarios", {
    
    # Test normal case with varied logFC values
    logFC_values <- c(-2, -1, 0, 1, 2)
    colors <- mapLogFCToColor(logFC_values)
    expect_equal(length(colors), 5)
    expect_true(all(grepl("^#[0-9A-Fa-f]{6}$", colors))) # Valid hex colors
    
    # Test case with all NA values
    na_values <- c(NA, NA, NA)
    na_colors <- mapLogFCToColor(na_values)
    expect_equal(length(na_colors), 3)
    expect_true(all(na_colors == "#D3D3D3"))
    
    # Test case with all same values
    same_values <- c(1, 1, 1)
    same_colors <- mapLogFCToColor(same_values)
    expect_equal(length(same_colors), 3)
    expect_true(all(same_colors == "#D3D3D3"))
    
    # Test empty input
    empty_colors <- mapLogFCToColor(numeric(0))
    expect_equal(length(empty_colors), 0)
})

# =============================================================================
# TESTS FOR RELATIONSHIP PROPERTIES
# =============================================================================

test_that("getRelationshipProperties returns correct structure", {
    props <- getRelationshipProperties()
    
    expect_type(props, "list")
    expect_true("complex" %in% names(props))
    expect_true("regulatory" %in% names(props))
    expect_true("phosphorylation" %in% names(props))
    expect_true("other" %in% names(props))
    
    # Test complex properties
    complex_props <- props$complex
    expect_true("types" %in% names(complex_props))
    expect_true("Complex" %in% complex_props$types)
    expect_equal(complex_props$consolidate, "undirected")
    
    # Test regulatory properties
    reg_props <- props$regulatory
    expect_true("colors" %in% names(reg_props))
    expect_true("Inhibition" %in% names(reg_props$colors))
    expect_equal(reg_props$consolidate, "bidirectional")
})

# =============================================================================
# TESTS FOR EDGE CONSOLIDATION
# =============================================================================

test_that("consolidateEdges properly consolidates bidirectional relationships", {
    
    # Create test edges with bidirectional regulatory relationships
    test_edges <- data.frame(
        source = c("TP53", "MDM2", "ATM", "BRCA1"),
        target = c("MDM2", "TP53", "TP53", "ATM"),
        interaction = c("Inhibition", "Inhibition", "Phosphorylation", "Complex"),
        stringsAsFactors = FALSE
    )
    
    consolidated <- consolidateEdges(test_edges)
    
    expect_s3_class(consolidated, "data.frame")
    expect_true("edge_type" %in% names(consolidated))
    expect_true("category" %in% names(consolidated))
    
    # Should have fewer edges than original due to consolidation
    expect_lt(nrow(consolidated), nrow(test_edges))
    
    # Check that bidirectional inhibition was consolidated
    inhibition_edges <- consolidated[grepl("Inhibition", consolidated$interaction), ]
    expect_equal(nrow(inhibition_edges), 1)
    expect_equal(inhibition_edges$edge_type, "bidirectional")
})

test_that("consolidateEdges handles empty input", {
    empty_edges <- data.frame(
        source = character(0),
        target = character(0),
        interaction = character(0),
        stringsAsFactors = FALSE
    )
    
    result <- consolidateEdges(empty_edges)
    expect_equal(nrow(result), 0)
})

# =============================================================================
# TESTS FOR EDGE STYLING
# =============================================================================

test_that("getEdgeStyle returns appropriate styling", {
    
    # Test regulatory relationship styling
    style <- getEdgeStyle("Inhibition", "regulatory", "directed")
    expect_type(style, "list")
    expect_true("color" %in% names(style))
    expect_equal(style$color, "#FF4444") # Red for inhibition
    
    # Test complex relationship styling
    complex_style <- getEdgeStyle("Complex", "complex", "undirected")
    expect_equal(complex_style$arrow, "none")
    expect_equal(complex_style$color, "#8B4513")
    
    # Test unknown relationship
    unknown_style <- getEdgeStyle("Unknown", "other", "directed")
    expect_equal(unknown_style$color, "#666666")
})

# =============================================================================
# TESTS FOR NODE ELEMENT CREATION
# =============================================================================

test_that("createNodeElements creates proper node structures", {
    nodes <- create_mock_subnetwork_nodes()
    
    # Test with default label type (id)
    node_elements <- createNodeElements(nodes, "id")
    expect_equal(length(node_elements), nrow(nodes))
    expect_true(all(grepl("data:", node_elements)))
    expect_true(all(grepl("id:", node_elements)))
    expect_true(all(grepl("label:", node_elements)))
    
    # Test with hgncName label type
    node_elements_hgnc <- createNodeElements(nodes, "hgncName")
    expect_equal(length(node_elements_hgnc), nrow(nodes))
    
    # Test nodes without logFC column
    nodes_no_logfc <- nodes[, !names(nodes) %in% "logFC"]
    node_elements_no_logfc <- createNodeElements(nodes_no_logfc, "id")
    expect_equal(length(node_elements_no_logfc), nrow(nodes_no_logfc))
})

# =============================================================================
# TESTS FOR EDGE ELEMENT CREATION
# =============================================================================

test_that("createEdgeElements creates proper edge structures", {
    edges <- create_mock_subnetwork_edges()
    
    edge_elements <- createEdgeElements(edges)
    expect_type(edge_elements, "list")
    expect_gt(length(edge_elements), 0)
    
    # Check that all elements contain required fields
    expect_true(all(sapply(edge_elements, function(x) grepl("source:", x))))
    expect_true(all(sapply(edge_elements, function(x) grepl("target:", x))))
    
    # Test empty edges
    empty_edges <- data.frame(
        source = character(0),
        target = character(0),
        interaction = character(0),
        stringsAsFactors = FALSE
    )
    empty_elements <- createEdgeElements(empty_edges)
    expect_equal(length(empty_elements), 0)
})

# =============================================================================
# TESTS FOR CYTOSCAPE CONFIG GENERATION
# =============================================================================

test_that("generateCytoscapeConfig creates complete configuration", {
    nodes <- create_mock_subnetwork_nodes()
    edges <- create_mock_subnetwork_edges()
    
    config <- generateCytoscapeConfig(nodes, edges)
    
    expect_type(config, "list")
    expect_true("elements" %in% names(config))
    expect_true("style" %in% names(config))
    expect_true("layout" %in% names(config))
    expect_true("container_id" %in% names(config))
    expect_true("js_code" %in% names(config))
    
    expect_equal(config$container_id, "network-cy")
    expect_type(config$js_code, "character")
    expect_gt(nchar(config$js_code), 100)
})

test_that("generateCytoscapeConfig accepts custom parameters", {
    nodes <- create_mock_subnetwork_nodes()
    edges <- create_mock_subnetwork_edges()
    
    custom_layout <- list(name = "grid", fit = FALSE)
    custom_handlers <- list(edge_click = "function() { console.log('test'); }")
    
    config <- generateCytoscapeConfig(
        nodes, 
        edges,
        container_id = "custom-container",
        event_handlers = custom_handlers,
        layout_options = custom_layout
    )
    
    expect_equal(config$container_id, "custom-container")
    expect_equal(config$layout$name, "grid")
    expect_false(config$layout$fit)
    expect_true(grepl("console.log", config$js_code))
})


# =============================================================================
# TESTS FOR STYLE CONVERSION FUNCTIONS
# =============================================================================

test_that("convertStyleToJS creates valid JavaScript", {
    style_list <- list(
        list(
            selector = "node",
            style = list(
                `background-color` = "data(color)",
                width = "60px"
            )
        )
    )
    
    js_style <- convertStyleToJS(style_list)
    expect_type(js_style, "character")
    expect_true(grepl("selector", js_style))
    expect_true(grepl("background-color", js_style))
    expect_true(grepl("data\\(color\\)", js_style))
})

test_that("convertLayoutToJS creates valid JavaScript", {
    layout_list <- list(
        name = "dagre",
        fit = TRUE,
        padding = 30
    )
    
    js_layout <- convertLayoutToJS(layout_list)
    expect_type(js_layout, "character")
    expect_true(grepl("\"name\": \"dagre\"", js_layout))
    expect_true(grepl("\"fit\": true", js_layout))
    expect_true(grepl("\"padding\": 30", js_layout))
})

test_that("createNodeElements handles different label types", {
    nodes <- create_mock_subnetwork_nodes()
    
    # Test with id labels
    elements_id <- createNodeElements(nodes, "id")
    expect_true(all(grepl("P53_HUMAN|MDM2_HUMAN|ATM_HUMAN|BRCA1_HUMAN", elements_id)))
    
    # Test with hgncName labels
    elements_hgnc <- createNodeElements(nodes, "hgncName")
    expect_true(all(grepl("TP53|MDM2|ATM|BRCA1", elements_hgnc)))
    
    # Test with nodes missing hgncName
    nodes_no_hgnc <- nodes
    nodes_no_hgnc$hgncName <- NA
    elements_fallback <- createNodeElements(nodes_no_hgnc, "hgncName")
    expect_true(all(grepl("P53_HUMAN|MDM2_HUMAN|ATM_HUMAN|BRCA1_HUMAN", elements_fallback)))
})