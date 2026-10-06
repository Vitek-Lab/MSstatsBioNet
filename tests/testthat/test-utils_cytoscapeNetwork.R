# =============================================================================
# MOCK DATA
# =============================================================================

create_mock_nodes <- function() {
    data.frame(
        id          = c("P53_HUMAN", "MDM2_HUMAN", "ATM_HUMAN", "BRCA1_HUMAN"),
        logFC       = c(2.5, -1.8, 1.2, -2.1),
        pvalue      = c(0.001, 0.02, 0.03, 0.005),
        entity_name = c("TP53", "MDM2", "ATM", "BRCA1"),
        stringsAsFactors = FALSE
    )
}

create_mock_nodes_ptm <- function() {
    data.frame(
        id          = c("P53_HUMAN", "MDM2_HUMAN"),
        logFC       = c(2.5, -1.8),
        entity_name = c("TP53", "MDM2"),
        site        = c(NA, "S15_S20"),
        stringsAsFactors = FALSE
    )
}

create_mock_edges <- function() {
    data.frame(
        source      = c("P53_HUMAN", "MDM2_HUMAN", "ATM_HUMAN", "P53_HUMAN", "BRCA1_HUMAN"),
        target      = c("MDM2_HUMAN", "P53_HUMAN",  "P53_HUMAN", "BRCA1_HUMAN", "P53_HUMAN"),
        interaction = c("Inhibition", "Inhibition", "Phosphorylation", "Complex", "Complex"),
        evidence_url = c("link1", "link2", "link3", "link4", "link5"),
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

test_that(".mapLogFCToColor returns the no-data colour for all-NA input", {
    colors <- MSstatsBioNet:::.mapLogFCToColor(c(NA, NA, NA))
    expect_length(colors, 3)
    expect_true(all(colors == MSstatsBioNet:::NO_DATA_NODE_COLOR))
})

test_that(".mapLogFCToColor gives NA a colour distinct from a zero fold change", {
    colors <- MSstatsBioNet:::.mapLogFCToColor(c(-2, 0, NA, 2))
    expect_equal(colors[3], MSstatsBioNet:::NO_DATA_NODE_COLOR)
    expect_equal(colors[2], MSstatsBioNet:::NEUTRAL_NODE_COLOR)
    expect_false(MSstatsBioNet:::NO_DATA_NODE_COLOR ==
                     MSstatsBioNet:::NEUTRAL_NODE_COLOR)
    # NA does not shift the colours of the other values
    expect_equal(colors[-3], MSstatsBioNet:::.mapLogFCToColor(c(-2, 0, 2)))
})

test_that(".mapLogFCToColor keeps NA distinct when one value is non-NA", {
    colors <- MSstatsBioNet:::.mapLogFCToColor(c(1, NA))
    expect_equal(colors, c(MSstatsBioNet:::NEUTRAL_NODE_COLOR,
                           MSstatsBioNet:::NO_DATA_NODE_COLOR))
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

test_that(".mapLogFCToColor handles Inf and -Inf values", {
    colors <- MSstatsBioNet:::.mapLogFCToColor(c(-Inf, 0, Inf))
    expect_length(colors, 3)
    expect_true(all(grepl("^#[0-9A-Fa-f]{6}$", colors)))
})

# =============================================================================
# .relProps
# =============================================================================

test_that(".relProps returns correct structure", {
    props <- MSstatsBioNet:::.relProps()
    
    expect_type(props, "list")
    expect_true(all(c("complex", "regulatory", "phosphorylation",
                      "modification", "other") %in% names(props)))
    
    expect_true("Inhibition" %in% names(props$regulatory$colors))
    expect_true("Activation" %in% names(props$regulatory$colors))
})

test_that(".relProps puts every contract statement type in exactly one category", {
    props <- MSstatsBioNet:::.relProps()
    listed <- unlist(lapply(props, `[[`, "types"), use.names = FALSE)
    expect_setequal(listed, MSstatsBioNet:::INTERACTION_TYPES)
    expect_false(anyDuplicated(listed) > 0)
})

# =============================================================================
# .classify
# =============================================================================

test_that(".classify maps interaction types to correct categories", {
    expect_equal(MSstatsBioNet:::.classify("Inhibition"),      "regulatory")
    expect_equal(MSstatsBioNet:::.classify("Activation"),      "regulatory")
    expect_equal(MSstatsBioNet:::.classify("Phosphorylation"), "phosphorylation")
    expect_equal(MSstatsBioNet:::.classify("Complex"),         "complex")
    expect_equal(MSstatsBioNet:::.classify("Association"),     "complex")
    expect_equal(MSstatsBioNet:::.classify("Dephosphorylation"), "phosphorylation")
    expect_equal(MSstatsBioNet:::.classify("Ubiquitination"),  "modification")
    expect_equal(MSstatsBioNet:::.classify("Translocation"),   "other")
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

test_that(".edgeStyle falls back to grey for unknown category", {
    style <- MSstatsBioNet:::.edgeStyle("Unknown", "other", "directed")
    expect_equal(style$color, "#666666")
})

# =============================================================================
# .consolidateEdges
# =============================================================================

test_that(".consolidateEdges keeps opposite regulatory edges as separate directed edges", {
    edges <- create_mock_edges()
    result <- MSstatsBioNet:::.consolidateEdges(edges)

    expect_s3_class(result, "data.frame")
    expect_true(all(c("edge_type", "category", "ptm_overlap") %in% names(result)))

    # Two Inhibition edges in opposite directions → both kept as directed
    inhibition <- result[grepl("Inhibition", result$interaction), ]
    expect_equal(nrow(inhibition), 2)
    expect_true(all(inhibition$edge_type == "directed"))
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

test_that(".consolidateEdges reads the directed column", {
    edges <- data.frame(source      = c("A", "B", "C"),
                        target      = c("B", "A", "D"),
                        interaction = c("Activation", "Activation", "Complex"),
                        directed    = c(FALSE, FALSE, TRUE),
                        stringsAsFactors = FALSE)
    result <- MSstatsBioNet:::.consolidateEdges(edges)
    # Undirected A-B pair is drawn once; a directed Complex keeps its arrow
    expect_equal(nrow(result), 2)
    expect_equal(result$edge_type[result$interaction == "Activation"],
                 "undirected")
    expect_equal(result$edge_type[result$interaction == "Complex"], "directed")
})

test_that(".consolidateEdges falls back to the statement type without directed", {
    edges <- data.frame(source      = c("A", "C", "E"),
                        target      = c("B", "D", "F"),
                        interaction = c("Complex", "Association", "Activation"),
                        stringsAsFactors = FALSE)
    result <- MSstatsBioNet:::.consolidateEdges(edges)
    # A Complex edge with no reverse edge is still undirected
    expect_equal(result$edge_type, c("undirected", "undirected", "directed"))

    edges$directed <- c(NA, FALSE, NA)
    expect_equal(MSstatsBioNet:::.consolidateEdges(edges)$edge_type,
                 c("undirected", "undirected", "directed"))
})

test_that(".consolidateEdges keeps undirected edges of different types", {
    edges <- data.frame(source      = c("A", "B"),
                        target      = c("B", "A"),
                        interaction = c("Complex", "Association"),
                        stringsAsFactors = FALSE)
    result <- MSstatsBioNet:::.consolidateEdges(edges)
    expect_equal(nrow(result), 2)
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

test_that(".buildElements uses entity_name label when requested", {
    nodes  <- create_mock_nodes()
    result <- MSstatsBioNet:::.buildElements(nodes, data.frame(), "entity_name")
    
    protein_nodes <- Filter(function(el) !is.null(el$data$node_type) &&
                                el$data$node_type == "protein", result)
    labels <- sapply(protein_nodes, function(el) el$data$label)
    expect_true(all(labels %in% c("TP53", "MDM2", "ATM", "BRCA1")))
})

test_that(".buildElements falls back to id when entity_name is NA", {
    nodes <- create_mock_nodes()
    nodes$entity_name <- NA
    result <- MSstatsBioNet:::.buildElements(nodes, data.frame(), "entity_name")
    
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

# Node data of the main (non-PTM, non-compound) nodes, keyed by id
.main_node_data <- function(elements) {
    main <- Filter(function(el) identical(el$data$node_type, "protein"),
                   elements)
    stats::setNames(lapply(main, `[[`, "data"),
                    vapply(main, function(el) el$data$id, character(1)))
}

# A network with (a) a significant protein, (b) a protein in the input with
# logFC near 0 that was not queried, and (c) a protein not in the input
.node_status_fixture <- function() {
    data.frame(
        id                = c("SIG", "FLAT", "LATENT"),
        entity_type       = "protein",
        measured          = c(TRUE, TRUE, FALSE),
        included_in_query = c(TRUE, FALSE, TRUE),
        logFC             = c(2.5, 0.01, NA),
        adj.pvalue        = c(0.001, 0.9, NA),
        stringsAsFactors  = FALSE
    )
}

test_that(".node_display_status distinguishes measured, not queried, and latent", {
    status <- MSstatsBioNet:::.node_display_status(.node_status_fixture())
    expect_equal(status, c("measured", "not_queried", "latent"))
})

test_that(".node_display_status treats missing status columns as measured", {
    nodes <- create_mock_nodes()
    expect_true(all(MSstatsBioNet:::.node_display_status(nodes) == "measured"))
    nodes$logFC[2] <- NA
    expect_equal(MSstatsBioNet:::.node_display_status(nodes)[2], "no_logfc")
    nodes$logFC <- NULL
    expect_true(all(MSstatsBioNet:::.node_display_status(nodes) == "measured"))
})

test_that(".buildElements styles the node-status fixture three ways", {
    nodes <- .main_node_data(MSstatsBioNet:::.buildElements(
        .node_status_fixture(), data.frame()))
    expect_equal(nodes$SIG$status, "measured")
    expect_equal(nodes$FLAT$status, "not_queried")
    expect_equal(nodes$LATENT$status, "latent")
    expect_equal(nodes$LATENT$color, MSstatsBioNet:::NO_DATA_NODE_COLOR)
    expect_false(nodes$FLAT$color == nodes$LATENT$color)
})

test_that(".buildElements colours a PTM protein from its protein-level row", {
    nodes <- data.frame(id    = c("P1", "P1", "P1"),
                        site  = c("S5", NA, "T9"),
                        logFC = c(-2, 1.5, 2),
                        stringsAsFactors = FALSE)
    elements <- MSstatsBioNet:::.buildElements(nodes, data.frame())
    protein <- .main_node_data(elements)$P1
    expected <- MSstatsBioNet:::.mapLogFCToColor(nodes$logFC)
    expect_equal(protein$color, expected[2])
    expect_equal(protein$status, "measured")
    ptm <- Filter(function(el) identical(el$data$node_type, "ptm"), elements)
    ptm_colors <- vapply(ptm, function(el) el$data$color, character(1))
    names(ptm_colors) <- vapply(ptm, function(el) el$data$label, character(1))
    expect_equal(unname(ptm_colors[c("S5", "T9")]), expected[c(1, 3)])
})

test_that(".buildElements draws a PTM protein without a protein-level row as a container", {
    elements <- MSstatsBioNet:::.buildElements(create_mock_nodes_ptm(),
                                               data.frame())
    mdm2 <- .main_node_data(elements)$MDM2_HUMAN
    expect_equal(mdm2$status, "sites_only")
    expect_equal(mdm2$color, MSstatsBioNet:::NO_DATA_NODE_COLOR)
    # The site nodes keep the colour of their row
    ptm <- Filter(function(el) identical(el$data$node_type, "ptm"), elements)
    site_color <- MSstatsBioNet:::.mapLogFCToColor(
        create_mock_nodes_ptm()$logFC)[2]
    expect_true(all(vapply(ptm, function(el) el$data$color, character(1)) ==
                        site_color))
    expect_true(all(vapply(ptm, function(el) el$data$status, character(1)) ==
                        "measured"))
})

test_that(".buildElements gives nodes a shape by entity_type", {
    nodes <- data.frame(id          = c("P", "M", "D", "F", "X", "S"),
                        entity_type = c("protein", "metabolite", "drug",
                                        "family", NA, "ptm_site"),
                        stringsAsFactors = FALSE)
    data <- .main_node_data(MSstatsBioNet:::.buildElements(nodes, data.frame()))
    shapes <- vapply(data, `[[`, character(1), "shape")
    expect_equal(unname(shapes[c("P", "M", "D", "F", "X", "S")]),
                 c("round-rectangle", "hexagon", "diamond", "barrel",
                   "round-rectangle", "round-rectangle"))
    expect_equal(data$X$entity_type, "")
    expect_equal(data$M$entity_type, "metabolite")
})

test_that(".buildElements gives every node the protein shape without entity_type", {
    data <- .main_node_data(MSstatsBioNet:::.buildElements(create_mock_nodes(),
                                                           data.frame()))
    expect_true(all(vapply(data, `[[`, character(1), "shape") ==
                        "round-rectangle"))
})

test_that("every ENTITY_TYPES value has a node shape", {
    expect_setequal(names(MSstatsBioNet:::NODE_SHAPES),
                    MSstatsBioNet:::ENTITY_TYPES)
})

test_that(".buildElements draws undirected edges without an arrow", {
    edges <- data.frame(source = c("P53_HUMAN", "MDM2_HUMAN"),
                        target = c("MDM2_HUMAN", "ATM_HUMAN"),
                        interaction = c("Association", "Activation"),
                        directed = c(FALSE, TRUE),
                        stringsAsFactors = FALSE)
    elements <- MSstatsBioNet:::.buildElements(create_mock_nodes(), edges)
    edge_data <- lapply(Filter(function(el) !is.null(el$data$source) &&
                                   !identical(el$data$edge_type,
                                              "ptm_attachment"),
                               elements), `[[`, "data")
    arrows <- vapply(edge_data, `[[`, character(1), "arrow_shape")
    names(arrows) <- vapply(edge_data, `[[`, character(1), "interaction")
    expect_equal(unname(arrows[c("Association", "Activation")]),
                 c("none", "triangle"))
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

# ----- displayLabelType = "entityName" deprecation -----

test_that("cytoscapeNetwork accepts the deprecated displayLabelType = 'entityName'", {
    nodes <- create_mock_nodes()
    expect_warning(
        widget <- cytoscapeNetwork(nodes, data.frame(), displayLabelType = "entityName"),
        "\"entityName\" is deprecated"
    )
    labels <- vapply(widget$x$elements, function(el) el$data$label, "")
    expect_equal(labels, nodes$entity_name)
})

test_that("cytoscapeNetwork labels nodes by entity_name without a warning", {
    nodes <- create_mock_nodes()
    expect_silent(
        widget <- cytoscapeNetwork(nodes, data.frame(), displayLabelType = "entity_name")
    )
    labels <- vapply(widget$x$elements, function(el) el$data$label, "")
    expect_equal(labels, nodes$entity_name)
})

test_that("cytoscapeNetwork rejects an unknown displayLabelType", {
    expect_error(
        cytoscapeNetwork(create_mock_nodes(), data.frame(), displayLabelType = "name"),
        "displayLabelType must be \"id\" or \"entity_name\""
    )
})

# ----- Edge contract columns (Phase 1b of the API refactor) -----

test_that(".buildElements passes evidence_url through to the widget edges", {
    nodes <- data.frame(id = c("A", "B"), stringsAsFactors = FALSE)
    edges <- data.frame(source = "A", target = "B", interaction = "Activation",
                        evidence_url = "https://db.indra.bio/statements/from_hash/1?format=html",
                        stringsAsFactors = FALSE)
    elements <- MSstatsBioNet:::.buildElements(nodes, edges)
    edge_data <- Filter(function(el) !is.null(el$data$source), elements)[[1]]$data
    expect_equal(edge_data$evidence_url,
                 "https://db.indra.bio/statements/from_hash/1?format=html")
    expect_null(edge_data$evidenceLink)
})
