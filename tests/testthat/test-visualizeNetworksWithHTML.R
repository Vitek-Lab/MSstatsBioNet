# =============================================================================
# MOCK DATA
# =============================================================================

create_mock_nodes <- function() {
    data.frame(
        id       = c("P53_HUMAN", "MDM2_HUMAN", "ATM_HUMAN", "BRCA1_HUMAN"),
        logFC    = c(2.5, -1.8, 1.2, -2.1),
        pvalue   = c(0.001, 0.02, 0.03, 0.005),
        hgncName = c("TP53", "MDM2", "ATM", "BRCA1"),
        stringsAsFactors = FALSE
    )
}

create_mock_nodes_ptm <- function() {
    data.frame(
        id       = c("P53_HUMAN", "MDM2_HUMAN"),
        logFC    = c(2.5, -1.8),
        hgncName = c("TP53", "MDM2"),
        Site     = c(NA, "S15_S20"),
        stringsAsFactors = FALSE
    )
}

create_mock_edges <- function() {
    data.frame(
        source      = c("P53_HUMAN", "MDM2_HUMAN", "ATM_HUMAN", "P53_HUMAN", "BRCA1_HUMAN"),
        target      = c("MDM2_HUMAN", "P53_HUMAN",  "P53_HUMAN", "BRCA1_HUMAN", "P53_HUMAN"),
        interaction = c("Inhibition", "Inhibition", "Phosphorylation", "Complex", "Complex"),
        evidenceLink = c("link1", "link2", "link3", "link4", "link5"),
        stringsAsFactors = FALSE
    )
}

create_mock_edges_ptm <- function() {
    data.frame(
        source      = c("P53_HUMAN"),
        target      = c("MDM2_HUMAN"),
        interaction = c("Phosphorylation"),
        site        = c("S15"),
        stringsAsFactors = FALSE
    )
}

# =============================================================================
# .mapLogFCToColor
# =============================================================================

test_that(".mapLogFCToColor returns valid hex colours", {
    colors <- MSstatsBioNet:::.mapLogFCToColor(c(-2, -1, 0, 1, 2))
    expect_length(colors, 5)
    expect_true(all(grepl("^#[0-9A-Fa-f]{6}$", colors)))
})

test_that(".mapLogFCToColor returns grey for all-NA input", {
    colors <- MSstatsBioNet:::.mapLogFCToColor(c(NA, NA, NA))
    expect_length(colors, 3)
    expect_true(all(colors == "#D3D3D3"))
})

test_that(".mapLogFCToColor returns grey for all-identical values", {
    colors <- MSstatsBioNet:::.mapLogFCToColor(c(1, 1, 1))
    expect_length(colors, 3)
    expect_true(all(colors == "#D3D3D3"))
})

test_that(".mapLogFCToColor handles empty input", {
    colors <- MSstatsBioNet:::.mapLogFCToColor(numeric(0))
    expect_length(colors, 0)
})

# =============================================================================
# .relProps
# =============================================================================

test_that(".relProps returns correct structure", {
    props <- MSstatsBioNet:::.relProps()
    
    expect_type(props, "list")
    expect_true(all(c("complex", "regulatory", "phosphorylation", "other") %in% names(props)))
    
    expect_equal(props$complex$consolidate, "undirected")
    expect_equal(props$regulatory$consolidate, "bidirectional")
    expect_equal(props$phosphorylation$consolidate, "directed")
    
    expect_true("Inhibition" %in% names(props$regulatory$colors))
    expect_true("Activation" %in% names(props$regulatory$colors))
})

# =============================================================================
# .classify
# =============================================================================

test_that(".classify maps interaction types to correct categories", {
    expect_equal(MSstatsBioNet:::.classify("Inhibition"),      "regulatory")
    expect_equal(MSstatsBioNet:::.classify("Activation"),      "regulatory")
    expect_equal(MSstatsBioNet:::.classify("Phosphorylation"), "phosphorylation")
    expect_equal(MSstatsBioNet:::.classify("Complex"),         "complex")
    expect_equal(MSstatsBioNet:::.classify("Unknown"),         "other")
})

# =============================================================================
# .edgeStyle
# =============================================================================

test_that(".edgeStyle returns correct colour for regulatory interactions", {
    style <- MSstatsBioNet:::.edgeStyle("Inhibition", "regulatory", "directed")
    expect_type(style, "list")
    expect_equal(style$color, "#FF4444")
    
    style_act <- MSstatsBioNet:::.edgeStyle("Activation", "regulatory", "directed")
    expect_equal(style_act$color, "#44AA44")
})

test_that(".edgeStyle returns no arrow for undirected complex edges", {
    style <- MSstatsBioNet:::.edgeStyle("Complex", "complex", "undirected")
    expect_equal(style$arrow, "none")
    expect_equal(style$color, "#8B4513")
})

test_that(".edgeStyle returns triangle arrows for bidirectional edges", {
    style <- MSstatsBioNet:::.edgeStyle("Inhibition (bidirectional)", "regulatory", "bidirectional")
    expect_equal(style$arrow, "triangle")
})

test_that(".edgeStyle falls back to grey for unknown category", {
    style <- MSstatsBioNet:::.edgeStyle("Unknown", "other", "directed")
    expect_equal(style$color, "#666666")
})

# =============================================================================
# .consolidateEdges
# =============================================================================

test_that(".consolidateEdges consolidates bidirectional inhibition into one edge", {
    edges <- create_mock_edges()
    result <- MSstatsBioNet:::.consolidateEdges(edges)
    
    expect_s3_class(result, "data.frame")
    expect_true(all(c("edge_type", "category", "ptm_overlap") %in% names(result)))
    
    # Two Inhibition edges in opposite directions → one bidirectional edge
    inhibition <- result[grepl("Inhibition", result$interaction), ]
    expect_equal(nrow(inhibition), 1)
    expect_equal(inhibition$edge_type, "bidirectional")
})

test_that(".consolidateEdges marks phosphorylation as directed", {
    edges <- create_mock_edges()
    result <- MSstatsBioNet:::.consolidateEdges(edges)
    
    phospho <- result[result$interaction == "Phosphorylation", ]
    expect_equal(nrow(phospho), 1)
    expect_equal(phospho$edge_type, "directed")
    expect_equal(phospho$category, "phosphorylation")
})

test_that(".consolidateEdges marks complex as undirected", {
    edges <- create_mock_edges()
    result <- MSstatsBioNet:::.consolidateEdges(edges)
    
    complex <- result[result$interaction == "Complex", ]
    expect_equal(nrow(complex), 1)
    expect_equal(complex$edge_type, "undirected")
})

test_that(".consolidateEdges handles empty input", {
    empty <- data.frame(source = character(0), target = character(0),
                        interaction = character(0), stringsAsFactors = FALSE)
    result <- MSstatsBioNet:::.consolidateEdges(empty)
    expect_equal(nrow(result), 0)
})

# =============================================================================
# .ptmOverlap
# =============================================================================

test_that(".ptmOverlap detects overlapping PTM sites", {
    nodes <- create_mock_nodes_ptm()
    edges <- create_mock_edges_ptm()
    
    result <- MSstatsBioNet:::.ptmOverlap(edges, nodes)
    expect_type(result, "character")
    expect_length(result, 1)
    expect_true(grepl("S15", result[[1]]))
})

test_that(".ptmOverlap returns empty string when no overlap", {
    nodes <- create_mock_nodes_ptm()
    edges <- data.frame(source = "P53_HUMAN", target = "MDM2_HUMAN",
                        interaction = "Phosphorylation", site = "T999",
                        stringsAsFactors = FALSE)
    result <- MSstatsBioNet:::.ptmOverlap(edges, nodes)
    expect_equal(result[[1]], "")
})

test_that(".ptmOverlap handles empty edges gracefully", {
    nodes  <- create_mock_nodes_ptm()
    empty  <- data.frame(source = character(0), target = character(0),
                         interaction = character(0), stringsAsFactors = FALSE)
    result <- MSstatsBioNet:::.ptmOverlap(empty, nodes)
    expect_length(result, 0)
})

# =============================================================================
# .buildElements
# =============================================================================

test_that(".buildElements returns a list of elements", {
    nodes  <- create_mock_nodes()
    edges  <- create_mock_edges()
    result <- MSstatsBioNet:::.buildElements(nodes, edges)
    
    expect_type(result, "list")
    expect_gt(length(result), nrow(nodes))  # nodes + edges
})

test_that(".buildElements assigns correct node_type to proteins", {
    nodes  <- create_mock_nodes()
    result <- MSstatsBioNet:::.buildElements(nodes, data.frame())
    
    node_types <- sapply(result, function(el) el$data$node_type)
    expect_true(all(node_types == "protein"))
})

test_that(".buildElements creates PTM child nodes and attachment edges", {
    nodes  <- create_mock_nodes_ptm()
    result <- MSstatsBioNet:::.buildElements(nodes, data.frame())
    
    node_types <- sapply(result, function(el) el$data$node_type)
    expect_true("ptm" %in% node_types)
    expect_true("compound" %in% node_types)
    
    edge_types <- sapply(result, function(el) el$data$edge_type)
    expect_true("ptm_attachment" %in% edge_types)
})

test_that(".buildElements uses hgncName label when requested", {
    nodes  <- create_mock_nodes()
    result <- MSstatsBioNet:::.buildElements(nodes, data.frame(), "hgncName")
    
    protein_nodes <- Filter(function(el) !is.null(el$data$node_type) &&
                                el$data$node_type == "protein", result)
    labels <- sapply(protein_nodes, function(el) el$data$label)
    expect_true(all(labels %in% c("TP53", "MDM2", "ATM", "BRCA1")))
})

test_that(".buildElements falls back to id when hgncName is NA", {
    nodes <- create_mock_nodes()
    nodes$hgncName <- NA
    result <- MSstatsBioNet:::.buildElements(nodes, data.frame(), "hgncName")
    
    protein_nodes <- Filter(function(el) !is.null(el$data$node_type) &&
                                el$data$node_type == "protein", result)
    labels <- sapply(protein_nodes, function(el) el$data$label)
    expect_true(all(labels %in% nodes$id))
})

test_that(".buildElements computes width and height from label length", {
    nodes  <- create_mock_nodes()
    result <- MSstatsBioNet:::.buildElements(nodes, data.frame())
    
    protein_nodes <- Filter(function(el) !is.null(el$data$node_type) &&
                                el$data$node_type == "protein", result)
    widths  <- sapply(protein_nodes, function(el) el$data$width)
    heights <- sapply(protein_nodes, function(el) el$data$height)
    
    expect_true(all(widths  >= 60  & widths  <= 150))
    expect_true(all(heights >= 40  & heights <= 60))
})

test_that(".buildElements uses grey when logFC column is absent", {
    nodes <- create_mock_nodes()[, !names(create_mock_nodes()) %in% "logFC"]
    result <- MSstatsBioNet:::.buildElements(nodes, data.frame())
    
    protein_nodes <- Filter(function(el) !is.null(el$data$node_type) &&
                                el$data$node_type == "protein", result)
    colors <- sapply(protein_nodes, function(el) el$data$color)
    expect_true(all(colors == "#D3D3D3"))
})

# =============================================================================
# cytoscapeNetwork() — public API
# =============================================================================

test_that("cytoscapeNetwork() returns an htmlwidget", {
    w <- cytoscapeNetwork(create_mock_nodes(), create_mock_edges())
    expect_s3_class(w, "htmlwidget")
    expect_s3_class(w, "cytoscapeNetwork")
})

test_that("cytoscapeNetwork() x list contains elements and layout", {
    w <- cytoscapeNetwork(create_mock_nodes(), create_mock_edges())
    expect_true("elements" %in% names(w$x))
    expect_true("layout"   %in% names(w$x))
    expect_gt(length(w$x$elements), 0)
})

test_that("cytoscapeNetwork() passes custom layout options through", {
    w <- cytoscapeNetwork(create_mock_nodes(), create_mock_edges(),
                          layoutOptions = list(rankDir = "LR", rankSep = 120))
    expect_equal(w$x$layout$rankDir, "LR")
    expect_equal(w$x$layout$rankSep, 120)
})

test_that("cytoscapeNetwork() passes nodeFontSize through", {
    w <- cytoscapeNetwork(create_mock_nodes(), create_mock_edges(), nodeFontSize = 18)
    expect_equal(w$x$node_font_size, 18)
})

test_that("cytoscapeNetwork() accepts empty edges", {
    empty_edges <- data.frame(source = character(0), target = character(0),
                              interaction = character(0), stringsAsFactors = FALSE)
    expect_no_error(cytoscapeNetwork(create_mock_nodes(), empty_edges))
})

test_that("cytoscapeNetwork() errors when nodes has no id column", {
    bad_nodes <- data.frame(name = c("A", "B"), stringsAsFactors = FALSE)
    expect_error(cytoscapeNetwork(bad_nodes), "`id` column")
})

test_that("cytoscapeNetwork() errors when nodes is not a data frame", {
    expect_error(cytoscapeNetwork(list(id = "A")), "id column|data frame")
})