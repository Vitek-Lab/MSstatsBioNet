# ----- merge_networks() (Phase 5a of the API refactor) -----

# A small network that meets the contract. Each argument is one value per
# edge; nodes are made for every source and target, measured, with logFC 1.
.make_merge_network <- function(source = "A", target = "B",
                                interaction = "Activation",
                                statement_id = "1", evidence_count = 3L,
                                confidence = 0.9,
                                backend_database = "INDRA",
                                query_type = "subnetwork",
                                node_role = "passed_cutoffs") {
    edges <- data.frame(
        source = source, target = target, interaction = interaction,
        directed = !interaction %in% UNDIRECTED_INTERACTION_TYPES,
        site = NA_character_, confidence = confidence,
        evidence_count = evidence_count,
        evidence_url = paste0("https://example.org/", statement_id),
        statement_id = statement_id, backend_database = backend_database,
        query_type = query_type, stringsAsFactors = FALSE
    )
    node_ids <- unique(c(source, target))
    nodes <- data.frame(
        id = node_ids, entity_type = "protein", entity_name = node_ids,
        namespace = "HGNC", entity_id = node_ids, measured = TRUE,
        included_in_query = TRUE, node_role = node_role,
        site = NA_character_, has_measured_sites = FALSE, logFC = 1,
        adj.pvalue = 0.01, stringsAsFactors = FALSE
    )
    list(nodes = nodes, edges = edges)
}

test_that("the same statement from two queries becomes one edge", {
    subnetwork <- .make_merge_network()
    mediated <- .make_merge_network(query_type = "mediated")
    network <- merge_networks(subnetwork, mediated)
    expect_equal(nrow(network$edges), 1)
    expect_equal(network$edges$query_type, "subnetwork;mediated")
    expect_equal(nrow(network$nodes), 2)
})

test_that("query_type values are not repeated when merged networks merge again", {
    merged <- merge_networks(.make_merge_network(),
                             .make_merge_network(query_type = "mediated"))
    again <- merge_networks(merged, .make_merge_network())
    expect_equal(again$edges$query_type, "subnetwork;mediated")
})

test_that("the same relation from two backends stays two edges", {
    indra <- .make_merge_network()
    other <- .make_merge_network(backend_database = "STRING", confidence = 0.5)
    network <- merge_networks(indra, other)
    expect_equal(nrow(network$edges), 2)
    expect_setequal(network$edges$backend_database, c("INDRA", "STRING"))
    expect_setequal(network$edges$confidence, c(0.9, 0.5))
})

test_that("an undirected statement listed in both directions keeps both rows", {
    complex <- .make_merge_network(source = c("A", "B"), target = c("B", "A"),
                                   interaction = "Complex",
                                   statement_id = c("7", "7"))
    network <- merge_networks(complex, complex)
    expect_equal(nrow(network$edges), 2)
    expect_equal(network$edges$source, c("A", "B"))
})

test_that("differing copies keep the one with the most evidence, with a message", {
    fewer <- .make_merge_network(evidence_count = 2L)
    more <- .make_merge_network(evidence_count = 5L, confidence = 0.95,
                                query_type = "mediated")
    expect_message(network <- merge_networks(fewer, more),
                   "differ in evidence_count or confidence, statement_id: 1")
    expect_equal(network$edges$evidence_count, 5L)
    expect_equal(network$edges$confidence, 0.95)
    expect_equal(network$edges$query_type, "subnetwork;mediated")
})

test_that("ties in evidence_count keep the first copy", {
    first <- .make_merge_network(confidence = 0.8)
    second <- .make_merge_network(confidence = 0.6)
    expect_message(network <- merge_networks(first, second), "statement_id")
    expect_equal(network$edges$confidence, 0.8)
})

test_that("identical copies give no message", {
    expect_no_message(merge_networks(.make_merge_network(),
                                     .make_merge_network()))
})

test_that("an NA evidence_count loses to a count", {
    no_count <- .make_merge_network()
    no_count$edges$evidence_count <- NA_integer_
    counted <- .make_merge_network(evidence_count = 4L)
    with_na <- suppressMessages(
        .merge_edges(list(no_count$edges, counted$edges)))
    expect_equal(with_na$evidence_count, 4L)
})

test_that("columns that only one network has are filled with NA", {
    with_extra <- .make_merge_network()
    with_extra$edges$interaction_raw <- "Activation"
    without <- .make_merge_network(source = "B", target = "C",
                                   statement_id = "2")
    network <- merge_networks(with_extra, without)
    expect_equal(network$edges$interaction_raw, c("Activation", NA))
    expect_type(network$edges$interaction_raw, "character")
})

test_that("node roles are joined and has_measured_sites is TRUE if in any network", {
    first <- .make_merge_network()
    second <- .make_merge_network(node_role = "mediator",
                                  query_type = "mediated")
    second$nodes$has_measured_sites[second$nodes$id == "B"] <- TRUE
    network <- merge_networks(first, second)
    expect_equal(network$nodes$node_role,
                 c("passed_cutoffs;mediator", "passed_cutoffs;mediator"))
    expect_equal(network$nodes$has_measured_sites, c(FALSE, TRUE))
})

test_that("PTM site rows of one protein stay separate nodes", {
    first <- .make_merge_network()
    first$nodes <- rbind(first$nodes, first$nodes[2, ])
    first$nodes$site <- c(NA, "S10", "T20")
    first$nodes$logFC <- c(1, 2, 3)
    second <- .make_merge_network(query_type = "mediated")
    second$nodes <- first$nodes
    network <- merge_networks(first, second)
    expect_equal(nrow(network$nodes), 3)
    expect_equal(network$nodes$site, c(NA, "S10", "T20"))
    expect_equal(network$nodes$logFC, c(1, 2, 3))
})

test_that("conflicting node status without entities is an error naming the node", {
    first <- .make_merge_network()
    second <- .make_merge_network(query_type = "mediated")
    second$nodes$logFC[second$nodes$id == "B"] <- 2
    expect_error(merge_networks(first, second),
                 "disagree on logFC for node\\(s\\) B.*entities =")

    second <- .make_merge_network(query_type = "mediated")
    second$nodes$measured[2] <- FALSE
    second$nodes$logFC[2] <- NA
    second$nodes$adj.pvalue[2] <- NA
    expect_error(merge_networks(first, second),
                 "disagree on measured, logFC, adj.pvalue for node\\(s\\) B,")
})

test_that("the error names the columns that differ across all conflicting nodes", {
    first <- .make_merge_network()
    second <- .make_merge_network(query_type = "mediated")
    second$nodes$logFC[second$nodes$id == "A"] <- 2
    second$nodes$adj.pvalue[second$nodes$id == "B"] <- 0.5
    expect_error(merge_networks(first, second),
                 "disagree on logFC, adj.pvalue for node\\(s\\) A, B,")
})

test_that("included_in_query may differ and is TRUE if the node was in any query", {
    subnetwork <- .make_merge_network()
    regulators <- .make_merge_network(node_role = "upstream_regulator",
                                      query_type = "upstream_regulators:shared")
    regulators$nodes$included_in_query[regulators$nodes$id == "A"] <- FALSE
    network <- merge_networks(regulators, subnetwork)
    expect_equal(network$nodes$included_in_query, c(TRUE, TRUE))

    only_latent <- merge_networks(regulators, regulators)
    expect_equal(only_latent$nodes$included_in_query, c(FALSE, TRUE))
})

test_that("entities recompute measured and statistics, and unmatched nodes become latent", {
    first <- .make_merge_network()
    second <- .make_merge_network(source = "B", target = "C",
                                  statement_id = "2")
    second$nodes$logFC <- 5
    second$nodes$included_in_query[second$nodes$id == "C"] <- FALSE
    second$nodes$measured[second$nodes$id == "C"] <- FALSE
    second$nodes$logFC[second$nodes$id == "C"] <- NA
    second$nodes$adj.pvalue[second$nodes$id == "C"] <- NA
    entities <- prepare_entities(
        data.frame(Protein = c("A", "B"), log2FC = c(-1, 2),
                   adj.pvalue = c(0.04, 0.001)),
        entity_type = "protein", id_type = "uniprot")
    expect_error(merge_networks(first, second), "entities =")
    network <- merge_networks(first, second, entities = entities)
    nodes <- network$nodes[order(network$nodes$id), ]
    expect_equal(nodes$id, c("A", "B", "C"))
    expect_equal(nodes$measured, c(TRUE, TRUE, FALSE))
    expect_equal(nodes$logFC, c(-1, 2, NA))
    expect_equal(nodes$adj.pvalue, c(0.04, 0.001, NA))
    expect_equal(nodes$included_in_query, c(TRUE, TRUE, FALSE))
})

test_that("entities match PTM site nodes by parent id and site", {
    network <- .make_merge_network()
    network$nodes$site[2] <- "S10"
    entities <- prepare_entities(
        data.frame(Protein = c("A", "B_S10"), log2FC = c(1, 3),
                   adj.pvalue = c(0.01, 0.02),
                   type = c("protein", "ptm_site")),
        entity_type = "type", id_type = "uniprot")
    merged <- merge_networks(network, entities = entities)
    expect_equal(merged$nodes$logFC, c(1, 3))
    expect_equal(merged$nodes$measured, c(TRUE, TRUE))
})

test_that("regulators and provenance tables are combined by row", {
    first <- .make_merge_network()
    first$regulators <- data.frame(id = "A", analysis = "discrete",
                                   regulator_q = 0.01)
    first$provenance <- data.frame(backend_database = "INDRA",
                                   query_type = "subnetwork")
    second <- .make_merge_network(query_type = "mediated")
    second$regulators <- data.frame(id = "X", analysis = "metabolite")
    second$provenance <- data.frame(backend_database = "INDRA",
                                    query_type = "mediated",
                                    backend_version = NA_character_)
    third <- .make_merge_network(query_type = "path")
    network <- merge_networks(first, second, third)
    expect_equal(network$regulators$analysis, c("discrete", "metabolite"))
    expect_equal(network$regulators$regulator_q, c(0.01, NA))
    expect_equal(nrow(network$provenance), 2)
    expect_equal(network$provenance$backend_version, c(NA_character_, NA))
})

test_that("networks without regulators or provenance give none", {
    network <- merge_networks(.make_merge_network(), .make_merge_network())
    expect_named(network, c("nodes", "edges"))
})

test_that("other network elements are dropped with a message", {
    first <- .make_merge_network()
    first$notes <- "made by hand"
    expect_message(network <- merge_networks(first), "dropped .* notes")
    expect_named(network, c("nodes", "edges"))
})

test_that("merge_networks() checks its input", {
    expect_error(merge_networks(), "at least one network")
    bad <- .make_merge_network()
    bad$edges$confidence <- 2
    expect_error(merge_networks(.make_merge_network(), bad),
                 "network 2: .*confidence")
    expect_error(merge_networks(good = .make_merge_network(), bad = bad),
                 "network 'bad'")
    expect_error(merge_networks(.make_merge_network(), entities = "A"),
                 "entities must be a data.frame")
})

test_that("validate_network() accepts ;-joined node roles and checks each", {
    network <- .make_merge_network(node_role = "passed_cutoffs;mediator")
    expect_no_error(validate_network(network))
    network$nodes$node_role[1] <- "passed_cutoffs;not_a_role"
    expect_error(validate_network(network), "not_a_role")
})

test_that("merging INDRA networks from one entity table gives their union", {
    local_mocked_bindings(
        .callIndraCogexApi = function(ns, ids, fio, cogex_url) {
            readRDS(system.file("extdata/indraResponse.rds",
                                package = "MSstatsBioNet"))
        })
    input <- data.table::fread(system.file("extdata/groupComparisonModel.csv",
                                           package = "MSstatsBioNet"))
    entities <- select_entities(.build_entities_from_annotated_input(input))
    indra <- indra_backend()
    full <- suppressMessages(get_network(indra, entities))
    complexes <- suppressMessages(
        get_network(indra, entities, interaction_types = "Complex"))
    others <- suppressMessages(
        get_network(indra, entities,
                    interaction_types = setdiff(INTERACTION_TYPES, "Complex")))
    network <- merge_networks(complexes, others)
    expect_equal(nrow(network$edges), nrow(full$edges))
    expect_setequal(paste(network$edges$source, network$edges$target,
                          network$edges$statement_id),
                    paste(full$edges$source, full$edges$target,
                          full$edges$statement_id))
    expect_setequal(network$nodes$id, full$nodes$id)
    expect_identical(merge_networks(full, complexes)$edges,
                     .reset_rownames(full$edges))
})
