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
        c("id", "entity_name", "namespace", "entity_id", "site", "logFC",
          "adj.pvalue")
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

# Phase 2 moves the INDRA code behind a backend object and must reproduce
# this output exactly. Regenerate the fixture with
# _fixtures/make_golden_subnetwork.R only when a change is intended.
test_that("getSubnetworkFromIndra reproduces the pinned golden output", {
    suppressWarnings(subnetwork <- .run_mocked_subnetwork(statement_types = NULL))
    golden <- readRDS(test_path("_fixtures", "golden_subnetwork.rds"))
    expect_identical(subnetwork, golden)
})
