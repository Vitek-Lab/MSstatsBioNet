# ----- filter_by_ptm_site() and add_entity_values() (Phase 5d) -----

# A PTM network: EGFR (P00533) with two measured sites, phosphorylated by
# SRC, a latent node with INDRA's name as its id, plus an isolated node
.make_ptm_network <- function(edge_sites = c("Y1092", "T5"),
                              edge_targets = rep("P00533",
                                                 length(edge_sites)),
                              interaction = "Phosphorylation") {
    n_edges <- length(edge_sites)
    list(
        nodes = data.frame(
            id = c("P00533", "P00533", "SRC", "LONE"),
            entity_type = c("ptm_site", "ptm_site", "protein", "protein"),
            entity_name = c("EGFR", "EGFR", "SRC", "LONE"),
            namespace = "HGNC", entity_id = c("3236", "3236", "11283", "9"),
            measured = c(TRUE, TRUE, FALSE, TRUE),
            included_in_query = TRUE,
            node_role = c("passed_cutoffs", "passed_cutoffs", "user_added",
                          "passed_cutoffs"),
            site = c("S1039_S1042", "Y1092", NA, NA),
            has_measured_sites = c(TRUE, TRUE, FALSE, FALSE),
            logFC = c(1.5, -0.4, NA, 2), adj.pvalue = c(0.01, 0.3, NA, 0.01),
            stringsAsFactors = FALSE),
        edges = data.frame(
            source = rep("SRC", n_edges), target = edge_targets,
            interaction = interaction, directed = TRUE, site = edge_sites,
            confidence = 0.9, evidence_count = 3L,
            evidence_url = paste0("https://example.org/", seq_len(n_edges)),
            statement_id = as.character(seq_len(n_edges)),
            backend_database = "INDRA", query_type = "subnetwork",
            stringsAsFactors = FALSE),
        provenance = data.frame(backend_database = "INDRA",
                                query_type = "subnetwork")
    )
}

# Global profiling results, grounded as convert_ids() would
.make_global_entities <- function(ids = c("P00533", "P12931"),
                                  namespace = c("HGNC", "HGNC"),
                                  entity_id = c("3236", "11283"),
                                  log2FC = seq_along(ids) / 10) {
    entities <- prepare_entities(
        data.frame(Protein = ids, log2FC = log2FC,
                   adj.pvalue = rep(0.05, length(ids))),
        entity_type = "protein", id_type = "uniprot")
    entities$namespace <- namespace
    entities$entity_id <- entity_id
    entities$entity_name <- ifelse(is.na(entity_id), NA_character_, ids)
    entities
}

test_that("filter_by_ptm_site() keeps edges whose site is a site of their target", {
    network <- .make_ptm_network(edge_sites = c("Y1092", "T5", "S1042"))
    expect_message(filtered <- filter_by_ptm_site(network),
                   "Dropping 1 edge")
    expect_equal(filtered$edges$site, c("Y1092", "S1042"))
    expect_silent(validate_network(filtered))
})

test_that("filter_by_ptm_site() checks each edge on its own, not per source-target pair", {
    # Same source, target, and interaction: the old pooled filter kept T5
    network <- .make_ptm_network(edge_sites = c("Y1092", "T5"))
    filtered <- suppressMessages(filter_by_ptm_site(network))
    expect_equal(filtered$edges$site, "Y1092")
    old <- .filterByPtmSite(network$nodes, network$edges, TRUE)
    expect_setequal(old$edges$site, c("Y1092", "T5"))
})

test_that("filter_by_ptm_site() splits node sites on ';' (protein groups)", {
    network <- .make_ptm_network(edge_sites = "T5")
    network$nodes$site[2] <- "Y1092;T5"
    filtered <- filter_by_ptm_site(network)
    expect_equal(filtered$edges$site, "T5")
})

test_that("filter_by_ptm_site() matches a site only on the edge's target", {
    network <- .make_ptm_network(edge_sites = c("Y1092", "Y1092"),
                                 edge_targets = c("P00533", "LONE"))
    filtered <- suppressMessages(filter_by_ptm_site(network))
    expect_equal(filtered$edges$target, "P00533")
})

test_that("filter_by_ptm_site() drops edges without a site", {
    network <- .make_ptm_network(edge_sites = c("Y1092", NA))
    network$edges$interaction[2] <- "Activation"
    filtered <- suppressMessages(filter_by_ptm_site(network))
    expect_equal(filtered$edges$interaction, "Phosphorylation")
})

test_that("filter_by_ptm_site() keeps isolated nodes and other elements", {
    network <- .make_ptm_network(edge_sites = c("Y1092", "T5"))
    filtered <- suppressMessages(filter_by_ptm_site(network))
    expect_true("LONE" %in% filtered$nodes$id)
    expect_identical(filtered$provenance, network$provenance)
})

test_that("filter_by_ptm_site() drops nodes left without edges", {
    network <- .make_ptm_network(edge_sites = c("Y1092", "Y1092"),
                                 edge_targets = c("P00533", "LONE"))
    network$nodes <- rbind(network$nodes, network$nodes[4, ])
    network$nodes$id[5] <- "OTHER"
    network$edges$source[2] <- "OTHER"
    filtered <- suppressMessages(filter_by_ptm_site(network))
    expect_false("OTHER" %in% filtered$nodes$id)
    # LONE lost its only edge too
    expect_false("LONE" %in% filtered$nodes$id)
})

test_that("filter_by_ptm_site() stops when no edges remain", {
    network <- .make_ptm_network(edge_sites = "T5")
    expect_error(filter_by_ptm_site(network), "No edges remain")
})

test_that("filter_by_ptm_site() stops on a network without PTM site rows", {
    network <- .make_ptm_network()
    network$nodes$site <- NA_character_
    network$nodes$has_measured_sites <- FALSE
    expect_error(filter_by_ptm_site(network), "no PTM site rows")
})

test_that("filter_by_ptm_site() checks the network contract", {
    network <- .make_ptm_network()
    network$edges$source <- NULL
    expect_error(filter_by_ptm_site(network), "contract")
})

test_that("add_entity_values() adds a protein-level row to a sites-only protein", {
    network <- .make_ptm_network(edge_sites = "Y1092")
    expect_message(result <- add_entity_values(network,
                                               .make_global_entities()),
                   "added a protein-level row to 1 protein")
    egfr <- result$nodes[result$nodes$id == "P00533", ]
    expect_equal(nrow(egfr), 3)
    protein_row <- egfr[is.na(egfr$site), ]
    expect_equal(protein_row$logFC, 0.1)
    expect_equal(protein_row$entity_type, "protein")
    expect_true(protein_row$measured)
    expect_true(protein_row$has_measured_sites)
    expect_equal(protein_row$node_role, "passed_cutoffs")
    # placed before the protein's site rows
    expect_equal(which(result$nodes$id == "P00533")[1],
                 which(result$nodes$id == "P00533" &
                           is.na(result$nodes$site)))
    # the site rows keep their own values
    expect_equal(egfr$logFC[!is.na(egfr$site)], c(1.5, -0.4))
    expect_silent(validate_network(result))
})

test_that("add_entity_values() matches a latent node by grounding", {
    network <- .make_ptm_network(edge_sites = "Y1092")
    result <- suppressMessages(add_entity_values(network,
                                                 .make_global_entities()))
    src <- result$nodes[result$nodes$id == "SRC", ]
    expect_equal(src$logFC, 0.2)
    expect_true(src$measured)
    expect_equal(src$node_role, "user_added")
})

test_that("add_entity_values() matches by grounding only rows not matched by id", {
    # Isoform node P04637-2 shares HGNC:11998 with P04637, which matches by id
    network <- .make_ptm_network(edge_sites = "Y1092")
    network$nodes$id[4] <- "P04637-2"
    network$nodes$entity_id[4] <- "11998"
    network$nodes[5, ] <- network$nodes[4, ]
    network$nodes$id[5] <- "P04637"
    network$nodes$logFC[4:5] <- NA
    network$nodes$adj.pvalue[4:5] <- NA
    entities <- .make_global_entities(ids = "P04637", namespace = "HGNC",
                                      entity_id = "11998")
    result <- suppressMessages(add_entity_values(network, entities))
    expect_equal(result$nodes$logFC[result$nodes$id == "P04637"], 0.1)
    expect_true(is.na(result$nodes$logFC[result$nodes$id == "P04637-2"]))
})

test_that("add_entity_values() gives no values to a node matching several rows", {
    network <- .make_ptm_network(edge_sites = "Y1092")
    entities <- .make_global_entities(ids = c("P12931", "P12931-2"),
                                      entity_id = c("11283", "11283"))
    suppressMessages(expect_message(result <- add_entity_values(network, entities),
                       "match several rows of entities.*SRC"))
    expect_true(is.na(result$nodes$logFC[result$nodes$id == "SRC"]))
    expect_false(result$nodes$measured[result$nodes$id == "SRC"])
})

test_that("add_entity_values() matches a node with several groundings on any of them", {
    network <- .make_ptm_network(edge_sites = "Y1092")
    network$nodes$namespace[3] <- "HGNC;STRING"
    network$nodes$entity_id[3] <- "11283;9606.ENSP1"
    entities <- .make_global_entities(ids = "P12931", namespace = "STRING",
                                      entity_id = "9606.ENSP1")
    result <- suppressMessages(add_entity_values(network, entities))
    expect_equal(result$nodes$logFC[result$nodes$id == "SRC"], 0.1)
})

test_that("add_entity_values() stops on nodes with values unless overwrite = TRUE", {
    network <- .make_ptm_network(edge_sites = "Y1092")
    entities <- .make_global_entities(ids = c("LONE", "P12931"),
                                      entity_id = c("9", "11283"))
    expect_error(add_entity_values(network, entities),
                 "LONE already have logFC")
    result <- suppressMessages(add_entity_values(network, entities,
                                                 overwrite = TRUE))
    expect_equal(result$nodes$logFC[result$nodes$id == "LONE"], 0.1)
})

test_that("add_entity_values() does not add entity rows that match no node", {
    network <- .make_ptm_network(edge_sites = "Y1092")
    entities <- .make_global_entities(ids = c("P12931", "Q99999"),
                                      entity_id = c("11283", "5"))
    result <- suppressMessages(add_entity_values(network, entities))
    expect_false("Q99999" %in% result$nodes$id)
    expect_equal(nrow(result$nodes), nrow(network$nodes))
})

test_that("add_entity_values() leaves edges and other elements unchanged", {
    network <- .make_ptm_network(edge_sites = "Y1092")
    result <- suppressMessages(add_entity_values(network,
                                                 .make_global_entities()))
    expect_identical(result$edges, network$edges)
    expect_identical(result$provenance, network$provenance)
})

test_that("add_entity_values() requires a grounded entity table", {
    network <- .make_ptm_network(edge_sites = "Y1092")
    entities <- .make_global_entities(namespace = c(NA, NA),
                                      entity_id = c(NA, NA))
    expect_error(add_entity_values(network, entities), "Run convert_ids()")
})

test_that("add_entity_values() reports ungrounded rows, which still match by id", {
    network <- .make_ptm_network(edge_sites = "Y1092")
    entities <- .make_global_entities(namespace = c(NA, "HGNC"),
                                      entity_id = c(NA, "11283"))
    suppressMessages(expect_message(result <- add_entity_values(network, entities),
                       "1 row\\(s\\) of entities have no grounding"))
    egfr <- result$nodes[result$nodes$id == "P00533" &
                             is.na(result$nodes$site), ]
    expect_equal(egfr$logFC, 0.1)
})

test_that("add_entity_values() copies site values from a PTM entity table", {
    network <- .make_ptm_network(edge_sites = "Y1092")
    network$nodes$logFC[2] <- NA
    network$nodes$adj.pvalue[2] <- NA
    entities <- prepare_entities(
        data.frame(Protein = "P00533_Y1092", log2FC = 3, adj.pvalue = 0.001),
        entity_type = "ptm_site", id_type = "uniprot")
    entities$namespace <- "HGNC"
    entities$entity_id <- "3236"
    entities$entity_name <- "EGFR"
    result <- suppressMessages(add_entity_values(network, entities))
    expect_equal(result$nodes$logFC[2], 3)
    # no protein-level entity row, so no protein-level node row is added
    expect_false(any(result$nodes$id == "P00533" & is.na(result$nodes$site)))
})

test_that("add_entity_values() checks its input", {
    network <- .make_ptm_network(edge_sites = "Y1092")
    entities <- .make_global_entities()
    expect_error(add_entity_values(network, entities, overwrite = NA),
                 "overwrite")
    expect_error(add_entity_values(network, data.frame(id = "A")),
                 "entities")
    bad <- network
    bad$nodes$id <- NULL
    expect_error(add_entity_values(bad, entities), "contract")
})

test_that("a data.table nodes table works in add_entity_values()", {
    network <- .make_ptm_network(edge_sites = "Y1092")
    network$nodes <- data.table::as.data.table(network$nodes)
    result <- suppressMessages(add_entity_values(network,
                                                 .make_global_entities()))
    expect_equal(nrow(result$nodes), 5)
})
