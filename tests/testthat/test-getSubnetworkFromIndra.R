test_that("getSubnetworkFromIndra works correctly", {
    input <- data.table::fread(
        system.file("extdata/groupComparisonModel.csv", package = "MSstatsBioNet")
    )
    local_mocked_bindings(.callIndraCogexApi = function(ns, ids, fio, cogex_url) {
        return(readRDS(system.file("extdata/indraResponse.rds", package = "MSstatsBioNet")))
    })
    suppressWarnings(subnetwork <- getSubnetworkFromIndra(input, statement_types = c("Activation", "Phosphorylation")))
    expect_equal(nrow(subnetwork$nodes), 3)
    expect_equal(nrow(subnetwork$edges), 2)
})

test_that("getSubnetworkFromIndra with different statement type works correctly", {
    input <- data.table::fread(
        system.file("extdata/groupComparisonModel.csv", package = "MSstatsBioNet")
    )
    local_mocked_bindings(.callIndraCogexApi = function(ns, ids, fio, cogex_url) {
        return(readRDS(system.file("extdata/indraResponse.rds", package = "MSstatsBioNet")))
    })
    suppressWarnings(
        subnetwork <- getSubnetworkFromIndra(input, statement_types = c("Complex"))
    )
    expect_equal(nrow(subnetwork$nodes), 8)
    expect_equal(nrow(subnetwork$edges), 16)
})

test_that("Exception is thrown for 400+ proteins in dataframe", {
    input_400 <- data.frame(
        Protein         = paste0("Protein", 1:400),
        log2FC          = rep(1.0, 400),
        adj.pvalue      = rep(0.05, 400),
        EntityNamespace = rep("HGNC", 400),
        EntityId        = paste0("HGNCID", 1:400),
        EntityName      = paste0("HGNCNAME", 1:400),
        issue           = NA
    )
    expect_error(
        getSubnetworkFromIndra(input_400),
        "Invalid Input Error: INDRA query must contain less than 400 proteins.  Consider lowering your p-value cutoff"
    )
})

test_that("Exception is thrown for missing columns in input", {
    input_missing_cols <- data.frame(
        Protein = paste0("Protein", 1:10),
        issue = NA,
        adj.pvalue = 0.05
    )
    expect_error(
        getSubnetworkFromIndra(input_missing_cols),
        "Invalid Input Error: input is missing required column\\(s\\): log2FC, EntityNamespace, EntityId, EntityName\\."
    )
})

# Build a subnetwork from the saved INDRA response, with the API call mocked
.run_mocked_subnetwork <- function(..., statement_types = c("Complex")) {
    input <- data.table::fread(
        system.file("extdata/groupComparisonModel.csv", package = "MSstatsBioNet")
    )
    local_mocked_bindings(
        .callIndraCogexApi = function(ns, ids, fio, cogex_url) {
            readRDS(system.file("extdata/indraResponse.rds", package = "MSstatsBioNet"))
        },
        .env = parent.frame()
    )
    getSubnetworkFromIndra(input, statement_types = statement_types, ...)
}

# ----- Deprecated arguments (Phase 0 of the API refactor) -----

test_that("paper_count_cutoff warns as deprecated and is ignored", {
    suppressWarnings(baseline <- .run_mocked_subnetwork())
    suppressWarnings(expect_warning(
        subnetwork <- .run_mocked_subnetwork(paper_count_cutoff = 1),
        "'paper_count_cutoff'.*deprecated"
    ))
    expect_equal(subnetwork, baseline)
    # Values >= 2 used to remove every edge and error; now they are ignored
    suppressWarnings(expect_warning(
        subnetwork <- .run_mocked_subnetwork(paper_count_cutoff = 2),
        "'paper_count_cutoff'.*deprecated"
    ))
    expect_equal(subnetwork, baseline)
})

test_that("correlation_cutoff warns as deprecated", {
    suppressWarnings(expect_warning(
        .run_mocked_subnetwork(correlation_cutoff = 0.5),
        "'correlation_cutoff'.*deprecated"
    ))
})

test_that("protein_level_data warns as deprecated and still adds correlations", {
    input <- data.table::fread(
        system.file("extdata/groupComparisonModel.csv", package = "MSstatsBioNet")
    )
    proteins <- unique(input$Protein)
    protein_level_data <- data.frame(
        Protein        = rep(proteins, each = 3),
        originalRUN    = rep(paste0("run", 1:3), times = length(proteins)),
        LogIntensities = sin(seq_len(3 * length(proteins)))
    )
    suppressWarnings(expect_warning(
        subnetwork <- .run_mocked_subnetwork(
            protein_level_data = protein_level_data, correlation_cutoff = 0
        ),
        "'protein_level_data'.*deprecated"
    ))
    expect_true("correlation" %in% colnames(subnetwork$edges))
})

test_that("default and explicit NULL arguments give no deprecation warning", {
    warns <- capture_warnings(.run_mocked_subnetwork())
    expect_false(any(grepl("deprecated", warns)))
    warns <- capture_warnings(.run_mocked_subnetwork(protein_level_data = NULL))
    expect_false(any(grepl("deprecated", warns)))
})

# ----- Edge and node contract (Phase 1 of the API refactor) -----

test_that("getSubnetworkFromIndra returns edges and nodes that meet the contract", {
    suppressWarnings(subnetwork <- .run_mocked_subnetwork())
    expect_silent(validate_network(subnetwork))
    expect_equal(
        colnames(subnetwork$edges),
        c("source", "target", "interaction", "directed", "site", "confidence",
          "evidence_count", "evidence_url", "statement_id", "backend_database",
          "query_type", "evidence_sources", "paperCount")
    )
    expect_equal(
        colnames(subnetwork$nodes),
        c("id", "entity_type", "entity_name", "namespace", "entity_id",
          "measured", "included_in_query", "node_role", "site",
          "has_measured_sites", "logFC", "adj.pvalue")
    )
})

test_that("getSubnetworkFromIndra maps INDRA fields onto the contract columns", {
    res <- readRDS(system.file("extdata/indraResponse.rds", package = "MSstatsBioNet"))
    suppressWarnings(subnetwork <- .run_mocked_subnetwork())
    edges <- subnetwork$edges

    expect_true(all(edges$interaction == "Complex"))
    expect_true(all(!edges$directed))
    expect_true(all(edges$backend_database == "INDRA"))
    expect_true(all(edges$query_type == "subnetwork"))
    expect_type(edges$evidence_count, "integer")
    expect_type(edges$statement_id, "character")

    # statement_id comes from stmt_json$matches_hash, which keeps full
    # precision, and confidence from data$belief
    by_hash <- lapply(res, function(stmt) {
        list(hash = jsonlite::fromJSON(stmt$data$stmt_json)$matches_hash,
             belief = stmt$data$belief)
    })
    hashes <- vapply(by_hash, function(x) x$hash, "")
    beliefs <- vapply(by_hash, function(x) x$belief, 1)
    expect_true(all(edges$statement_id %in% hashes))
    expect_equal(edges$confidence, unname(beliefs[match(edges$statement_id, hashes)]))

    expect_equal(
        edges$evidence_url,
        paste0("https://db.indra.bio/statements/from_hash/",
               edges$statement_id, "?format=html")
    )
})

test_that("getSubnetworkFromIndra marks only symmetric statement types undirected", {
    suppressWarnings(subnetwork <- .run_mocked_subnetwork(
        statement_types = c("Activation", "IncreaseAmount", "DecreaseAmount")
    ))
    expect_true(all(subnetwork$edges$directed))
})

test_that("getSubnetworkFromIndra returns character node IDs when fread reads them as integers", {
    suppressWarnings(subnetwork <- .run_mocked_subnetwork())
    expect_type(subnetwork$nodes$entity_id, "character")
})

test_that("getSubnetworkFromIndra subtracts incorrect curations when filter_by_curation = TRUE", {
    local_mocked_bindings(.get_incorrect_curation_count = function(stmt_hash) 1)
    suppressWarnings({
        uncurated <- .run_mocked_subnetwork(statement_types = "Activation")
        curated <- .run_mocked_subnetwork(statement_types = "Activation",
                                          filter_by_curation = TRUE)
    })
    kept <- uncurated$edges$evidence_count - 1L >= 1L
    expect_equal(curated$edges$statement_id, uncurated$edges$statement_id[kept])
    expect_equal(curated$edges$evidence_count,
                 uncurated$edges$evidence_count[kept] - 1L)
    expect_true(all(curated$nodes$id %in%
                    c(curated$edges$source, curated$edges$target)))
})

# ----- Golden output (pinned after Phase 1 of the API refactor) -----

# Later phases must reproduce this output exactly. Regenerate the fixture
# with _fixtures/make_golden_subnetwork.R only when a change is intended
# (last: Phase 3c, node status columns).
test_that("getSubnetworkFromIndra reproduces the pinned golden output", {
    suppressWarnings(subnetwork <- .run_mocked_subnetwork(statement_types = NULL))
    golden <- readRDS(test_path("_fixtures", "golden_subnetwork.rds"))
    expect_identical(subnetwork, golden)
})

# ----- Node status (Phase 3c of the API refactor) -----

# A minimal INDRA statement, as .callIndraCogexApi() returns it
.make_statement <- function(source, target, stmt_type = "Activation",
                            hash = "1", evidence_count = 1L) {
    list(
        data = list(belief = 0.9, evidence_count = evidence_count,
                    source_counts = '{"reach": 1}',
                    stmt_json = paste0('{"matches_hash": "', hash, '"}'),
                    stmt_type = stmt_type),
        source_ns = source[1], source_id = source[2], source_name = source[3],
        target_ns = target[1], target_id = target[2], target_name = target[3]
    )
}

# Annotated input: Protein, statistics, and groundings
.make_annotated_input <- function(Protein, EntityNamespace, EntityId,
                                  log2FC = rep(1, length(Protein)),
                                  adj.pvalue = rep(0.01, length(Protein)),
                                  ...) {
    data.frame(Protein = Protein, log2FC = log2FC, adj.pvalue = adj.pvalue,
               EntityNamespace = EntityNamespace, EntityId = EntityId,
               EntityName = paste0("name_", Protein), ...,
               stringsAsFactors = FALSE)
}

.run_with_statements <- function(input, statements, ...) {
    local_mocked_bindings(
        .callIndraCogexApi = function(ns, ids, fio, cogex_url) statements,
        .env = parent.frame()
    )
    suppressWarnings(getSubnetworkFromIndra(input, ...))
}

test_that("nodes in the input are measured, with their statistics", {
    suppressWarnings(subnetwork <- .run_mocked_subnetwork())
    nodes <- subnetwork$nodes
    input <- data.table::fread(
        system.file("extdata/groupComparisonModel.csv", package = "MSstatsBioNet")
    )
    expect_true(all(nodes$measured))
    expect_true(all(nodes$included_in_query))
    expect_true(all(nodes$node_role == "passed_cutoffs"))
    expect_true(all(nodes$entity_type == "protein"))
    expect_false(any(nodes$has_measured_sites))
    expect_equal(nodes$logFC, input$log2FC[match(nodes$id, input$Protein)])
    expect_equal(nodes$adj.pvalue,
                 input$adj.pvalue[match(nodes$id, input$Protein)])
})

test_that("force_include_other outside the input gives latent nodes with NA statistics", {
    input <- .make_annotated_input(c("A", "B"), "HGNC", c("1", "2"))
    statements <- list(
        .make_statement(c("HGNC", "1", "GENEA"), c("HGNC", "2", "GENEB"), hash = "1"),
        .make_statement(c("HGNC", "1", "GENEA"), c("HGNC", "99", "GENEZ"), hash = "2"),
        .make_statement(c("FPLX", "AKT", "AKT"), c("HGNC", "2", "GENEB"), hash = "3")
    )
    subnetwork <- .run_with_statements(input, statements,
                                       force_include_other = c("HGNC:99", "FPLX:AKT"))
    nodes <- subnetwork$nodes
    expect_setequal(nodes$id, c("A", "B", "GENEZ", "AKT"))
    latent <- nodes[nodes$id %in% c("GENEZ", "AKT"), ]
    expect_false(any(latent$measured))
    expect_true(all(is.na(latent$logFC)))
    expect_true(all(is.na(latent$adj.pvalue)))
    expect_true(all(latent$included_in_query))
    expect_true(all(latent$node_role == "user_added"))
    expect_equal(latent$entity_type[latent$id == "AKT"], "family")
    expect_equal(latent$entity_type[latent$id == "GENEZ"], "protein")
    expect_equal(latent$namespace[latent$id == "GENEZ"], "HGNC")
    expect_equal(latent$entity_id[latent$id == "GENEZ"], "99")
    expect_equal(latent$entity_name[latent$id == "GENEZ"], "GENEZ")
})

test_that("force_include_other in the input keeps the row's statistics, as user_added", {
    input <- .make_annotated_input(c("A", "B"), "HGNC", c("1", "2"),
                                   log2FC = c(2, 0.01),
                                   adj.pvalue = c(0.001, 0.9))
    statements <- list(
        .make_statement(c("HGNC", "1", "GENEA"), c("HGNC", "2", "GENEB")))
    subnetwork <- .run_with_statements(input, statements, pvalueCutoff = 0.05,
                                       force_include_other = "HGNC:2")
    nodes <- subnetwork$nodes
    b <- nodes[nodes$id == "B", ]
    expect_true(b$measured)
    expect_equal(b$logFC, 0.01)
    expect_equal(b$adj.pvalue, 0.9)
    expect_equal(b$node_role, "user_added")
    expect_equal(nodes$node_role[nodes$id == "A"], "passed_cutoffs")
})

test_that("a node is matched to the rows in the query before other rows", {
    # A2 grounds like A but fails the cutoff, so it doesn't make A ambiguous
    input <- .make_annotated_input(c("A", "A2", "B"), "HGNC", c("1", "1", "2"),
                                   adj.pvalue = c(0.01, 0.9, 0.01))
    statements <- list(
        .make_statement(c("HGNC", "1", "GENEA"), c("HGNC", "2", "GENEB")))
    subnetwork <- .run_with_statements(input, statements, pvalueCutoff = 0.05)
    expect_equal(subnetwork$edges$source, "A")
    expect_setequal(subnetwork$nodes$id, c("A", "B"))
})

test_that("a node matching rows of several nodes keeps INDRA's name and NA statistics", {
    input <- .make_annotated_input(c("A1", "A2", "B"), "HGNC", c("1", "1", "2"))
    statements <- list(
        .make_statement(c("HGNC", "1", "GENEA"), c("HGNC", "2", "GENEB")))
    expect_message(
        subnetwork <- .run_with_statements(input, statements),
        "1 node\\(s\\) from the backend match entities of several nodes.*GENEA")
    genea <- subnetwork$nodes[subnetwork$nodes$id == "GENEA", ]
    expect_true(genea$measured)
    expect_true(is.na(genea$logFC))
    expect_equal(genea$node_role, "passed_cutoffs")
    expect_silent(validate_network(subnetwork))
})

test_that("statements between identifiers shared by two namespaces stay separate edges", {
    # HGNC:2 and CHEBI:2 are different entities with the same identifier
    input <- .make_annotated_input(c("A", "B", "glucose"),
                                   c("HGNC", "HGNC", "CHEBI"), c("1", "2", "2"))
    statements <- list(
        .make_statement(c("HGNC", "1", "GENEA"), c("HGNC", "2", "GENEB"), hash = "1"),
        .make_statement(c("HGNC", "1", "GENEA"), c("CHEBI", "2", "glucose"), hash = "2")
    )
    subnetwork <- .run_with_statements(input, statements)
    expect_equal(nrow(subnetwork$edges), 2)
    expect_setequal(subnetwork$edges$target, c("B", "glucose"))
    expect_equal(subnetwork$nodes$entity_type[subnetwork$nodes$id == "glucose"],
                 "metabolite")
})

test_that("PTM site rows become rows of the protein's node, with has_measured_sites", {
    input <- .make_annotated_input(c("P1_S10", "P1_S20", "P2"), "HGNC",
                                   c("1", "1", "2"), log2FC = c(1, 2, 3),
                                   GlobalProtein = c("P1", "P1", "P2"))
    statements <- list(
        .make_statement(c("HGNC", "2", "GENEB"), c("HGNC", "1", "GENEA")))
    subnetwork <- .run_with_statements(input, statements)
    nodes <- subnetwork$nodes
    p1 <- nodes[nodes$id == "P1", ]
    expect_equal(p1$site, c("S10", "S20"))
    expect_equal(p1$logFC, c(1, 2))
    expect_true(all(p1$entity_type == "ptm_site"))
    expect_true(all(p1$has_measured_sites))
    expect_true(all(p1$measured))
    expect_false(nodes$has_measured_sites[nodes$id == "P2"])
    expect_equal(subnetwork$edges$target, "P1")
})

test_that("rows in the query with no grounding are dropped with a message", {
    input <- .make_annotated_input(c("A", "B", "C"), c("HGNC", "HGNC", NA),
                                   c("1", "2", NA))
    statements <- list(
        .make_statement(c("HGNC", "1", "GENEA"), c("HGNC", "2", "GENEB")))
    expect_message(.run_with_statements(input, statements),
                   "Dropping 1 row\\(s\\) with no entity grounding")
})

test_that("getSubnetworkFromIndra stops on input with several comparisons", {
    input <- .make_annotated_input(c("A", "B"), "HGNC", c("1", "2"),
                                   Label = c("T vs C", "U vs C"))
    expect_error(getSubnetworkFromIndra(input),
                 "input has 2 comparisons in its Label column: T vs C, U vs C")
})

test_that("getSubnetworkFromIndra checks force_include_other's type", {
    input <- .make_annotated_input(c("A", "B"), "HGNC", c("1", "2"))
    expect_error(getSubnetworkFromIndra(input, force_include_other = 1),
                 "force_include_other must be a character vector")
})
