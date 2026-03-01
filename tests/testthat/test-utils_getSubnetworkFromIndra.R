make_nodes <- function() {
    data.frame(
        id       = c("P53_HUMAN", "MDM2_HUMAN", "ATM_HUMAN"),
        logFC    = c(1.5, -1.0, 0.5),
        Site     = c("S15_S20", "T68",  NA),
        stringsAsFactors = FALSE
    )
}

make_edges <- function() {
    data.frame(
        source      = c("ATM_HUMAN",  "ATM_HUMAN",  "P53_HUMAN",  "P53_HUMAN"),
        target      = c("P53_HUMAN",  "MDM2_HUMAN", "MDM2_HUMAN", "ATM_HUMAN"),
        interaction = c("Phosphorylation", "Phosphorylation", "Activation", "Activation"),
        site        = c("S15",         "T999",        NA,            "S15"),
        stringsAsFactors = FALSE
    )
}

test_that(".filterByPtmSite returns input unchanged when filter_by_ptm_site = FALSE", {
    nodes  <- make_nodes()
    edges  <- make_edges()
    result <- MSstatsBioNet:::.filterByPtmSite(nodes, edges, filter_by_ptm_site = FALSE)
    
    expect_equal(nrow(result$nodes), nrow(nodes))
    expect_equal(nrow(result$edges), nrow(edges))
    expect_equal(result$nodes, nodes)
    expect_equal(result$edges, edges)
})

test_that(".filterByPtmSite returns input unchanged when no nodes have Site data", {
    nodes       <- make_nodes()
    nodes$Site  <- NA   # wipe all sites
    edges       <- make_edges()
    
    result <- MSstatsBioNet:::.filterByPtmSite(nodes, edges, filter_by_ptm_site = TRUE)
    
    expect_equal(nrow(result$nodes), nrow(nodes))
    expect_equal(nrow(result$edges), nrow(edges))
})

test_that(".filterByPtmSite keeps only edges with PTM site overlap on target node", {
    result <- MSstatsBioNet:::.filterByPtmSite(make_nodes(), make_edges(),
                                                  filter_by_ptm_site = TRUE)
    # Only ATM→P53 with site=S15 should survive
    expect_equal(nrow(result$edges), 1)
    expect_equal(result$edges$source,      "ATM_HUMAN")
    expect_equal(result$edges$target,      "P53_HUMAN")
    expect_equal(result$edges$site,        "S15")
})

test_that(".filterByPtmSite drops edges where edge site does not overlap node Site", {
    result <- MSstatsBioNet:::.filterByPtmSite(make_nodes(), make_edges(),
                                                  filter_by_ptm_site = TRUE)
    # ATM→MDM2 with site=T999 should be gone (T999 not in MDM2's T68)
    dropped <- result$edges[result$edges$source == "ATM_HUMAN" &
                                result$edges$target == "MDM2_HUMAN", ]
    expect_equal(nrow(dropped), 0)
})

test_that(".filterByPtmSite drops edges with NA edge site even if node has Site data", {
    result <- MSstatsBioNet:::.filterByPtmSite(make_nodes(), make_edges(),
                                                  filter_by_ptm_site = TRUE)
    # P53→MDM2 has site=NA so must be dropped
    dropped <- result$edges[result$edges$source == "P53_HUMAN" &
                                result$edges$target == "MDM2_HUMAN", ]
    expect_equal(nrow(dropped), 0)
})

test_that(".filterByPtmSite drops edges where target node has no Site data", {
    result <- MSstatsBioNet:::.filterByPtmSite(make_nodes(), make_edges(),
                                                  filter_by_ptm_site = TRUE)
    # P53→ATM: ATM node Site is NA so no overlap possible
    dropped <- result$edges[result$edges$source == "P53_HUMAN" &
                                result$edges$target == "ATM_HUMAN", ]
    expect_equal(nrow(dropped), 0)
})

test_that(".filterByPtmSite prunes nodes to only those in surviving edges", {
    result <- MSstatsBioNet:::.filterByPtmSite(make_nodes(), make_edges(),
                                                  filter_by_ptm_site = TRUE)
    # Only ATM_HUMAN (source) and P53_HUMAN (target) should remain
    expect_setequal(result$nodes$id, c("ATM_HUMAN", "P53_HUMAN"))
    expect_false("MDM2_HUMAN" %in% result$nodes$id)
})

test_that(".filterByPtmSite preserves all node columns after pruning", {
    result <- MSstatsBioNet:::.filterByPtmSite(make_nodes(), make_edges(),
                                                  filter_by_ptm_site = TRUE)
    expect_true(all(c("id", "logFC", "Site") %in% names(result$nodes)))
})

test_that(".filterByPtmSite keeps edge when site matches any of multiple node sites", {
    nodes <- data.frame(
        id   = c("A", "B"),
        Site = c("S15_S20_T68", NA),
        stringsAsFactors = FALSE
    )
    edges <- data.frame(
        source      = "A",
        target      = "B",
        interaction = "Phosphorylation",
        site        = "S20",    # matches second site in A's Site string
        stringsAsFactors = FALSE
    )
    # Note: filter checks target node — B has no Site, so this should drop.
    # Swap so A is the target to test multi-site matching.
    edges2 <- data.frame(
        source      = "B",
        target      = "A",
        interaction = "Phosphorylation",
        site        = "S20",
        stringsAsFactors = FALSE
    )
    result <- MSstatsBioNet:::.filterByPtmSite(nodes, edges2, filter_by_ptm_site = TRUE)
    expect_equal(nrow(result$edges), 1)
    expect_equal(result$edges$site, "S20")
})

test_that(".filterByPtmSite handles empty edges gracefully", {
    nodes <- make_nodes()
    empty_edges <- data.frame(
        source = character(0), target = character(0),
        interaction = character(0), site = character(0),
        stringsAsFactors = FALSE
    )
    result <- MSstatsBioNet:::.filterByPtmSite(nodes, empty_edges,
                                                  filter_by_ptm_site = TRUE)
    expect_equal(nrow(result$edges), 0)
    # All nodes pruned since no edges reference them
    expect_equal(nrow(result$nodes), 0)
})

test_that(".filterByPtmSite handles empty nodes gracefully", {
    empty_nodes <- data.frame(
        id = character(0), Site = character(0),
        stringsAsFactors = FALSE
    )
    edges <- make_edges()
    # No nodes have Site data so passthrough expected
    result <- MSstatsBioNet:::.filterByPtmSite(empty_nodes, edges,
                                                  filter_by_ptm_site = TRUE)
    expect_equal(nrow(result$edges), nrow(edges))
})

test_that(".filterByPtmSite always returns a list with nodes and edges", {
    result <- MSstatsBioNet:::.filterByPtmSite(make_nodes(), make_edges(),
                                                  filter_by_ptm_site = TRUE)
    expect_type(result, "list")
    expect_true(all(c("nodes", "edges") %in% names(result)))
    expect_s3_class(result$nodes, "data.frame")
    expect_s3_class(result$edges, "data.frame")
})