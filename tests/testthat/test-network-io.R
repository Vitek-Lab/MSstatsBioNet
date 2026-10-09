# ----- Provenance, save_network(), load_network() (Phase 5c) -----

.retrieval_time <- as.POSIXct("2026-10-09 12:00:00", tz = "UTC")

.mock_indra_with_clock <- function(env = parent.frame()) {
    local_mocked_bindings(
        .callIndraCogexApi = function(ns, ids, fio, cogex_url) {
            readRDS(system.file("extdata/indraResponse.rds",
                                package = "MSstatsBioNet"))
        },
        .get_current_time = function() .retrieval_time,
        .env = env
    )
}

.indra_input <- function() {
    input <- data.table::fread(system.file("extdata/groupComparisonModel.csv",
                                           package = "MSstatsBioNet"))
    select_entities(.build_entities_from_annotated_input(input))
}

.get_indra_network <- function(...) {
    suppressMessages(get_network(indra_backend(), .indra_input(), ...))
}

test_that("get_network() records one provenance row for INDRA", {
    .mock_indra_with_clock()
    network <- .get_indra_network(interaction_types = "Complex",
                                  min_evidence = 2)
    provenance <- network$provenance
    expect_s3_class(provenance, "data.frame")
    expect_named(provenance, PROVENANCE_COLUMNS)
    expect_equal(nrow(provenance), 1)
    expect_equal(provenance$backend_database, "INDRA")
    expect_equal(provenance$query_type, "subnetwork")
    expect_identical(provenance$retrieved_at, .retrieval_time)
    expect_identical(provenance$backend_version, NA_character_)
    expect_equal(provenance$backend_url, "https://discovery.indra.bio")
    expect_equal(provenance$organism, "9606")
    expect_equal(provenance$package_version,
                 as.character(utils::packageVersion("MSstatsBioNet")))
    parameters <- jsonlite::fromJSON(provenance$parameters)
    expect_equal(parameters$interaction_types, "Complex")
    expect_equal(parameters$min_evidence, 2)
    expect_null(parameters$evidence_sources)
    expect_true("evidence_sources" %in% names(parameters))
})

test_that("provenance joins to every edge on backend_database and query_type", {
    .mock_indra_with_clock()
    network <- .get_indra_network()
    expect_true(all(
        paste(network$edges$backend_database, network$edges$query_type) %in%
            paste(network$provenance$backend_database,
                  network$provenance$query_type)))
})

test_that("retrieved_at is taken when the response arrives", {
    calls <- character(0)
    local_mocked_bindings(
        .callIndraCogexApi = function(ns, ids, fio, cogex_url) {
            calls <<- c(calls, "request")
            readRDS(system.file("extdata/indraResponse.rds",
                                package = "MSstatsBioNet"))
        },
        .get_current_time = function() {
            calls <<- c(calls, "clock")
            .retrieval_time
        })
    .get_indra_network()
    expect_equal(calls[seq_len(2)], c("request", "clock"))
})

test_that("the backend's cogex_url is recorded", {
    .mock_indra_with_clock()
    network <- suppressMessages(get_network(
        indra_backend(cogex_url = "https://cogex.example.org"),
        .indra_input()))
    expect_equal(network$provenance$backend_url, "https://cogex.example.org")
})

test_that(".get_current_time() gives a UTC time", {
    now <- .get_current_time()
    expect_s3_class(now, "POSIXct")
    expect_equal(attr(now, "tzone"), "UTC")
    expect_lt(abs(as.numeric(difftime(now, Sys.time(), units = "secs"))), 5)
})

test_that("merge_networks() keeps the provenance of each network", {
    .mock_indra_with_clock()
    first <- .get_indra_network(interaction_types = "Complex")
    second <- .get_indra_network()
    second$provenance$retrieved_at <- .retrieval_time + 3600
    network <- merge_networks(first, second)
    expect_equal(nrow(network$provenance), 2)
    expect_s3_class(network$provenance$retrieved_at, "POSIXct")
    expect_identical(network$provenance$retrieved_at,
                     c(.retrieval_time, .retrieval_time + 3600))
})

test_that("filter_by_curation() keeps provenance", {
    .mock_indra_with_clock()
    network <- .get_indra_network()
    local_mocked_bindings(
        .count_incorrect_evidence = function(edges, backend) {
            integer(nrow(edges))
        })
    filtered <- filter_by_curation(network)
    expect_identical(filtered$provenance, network$provenance)
})

test_that("getSubnetworkFromIndra() output has no provenance, as before", {
    .mock_indra_with_clock()
    input <- data.table::fread(system.file("extdata/groupComparisonModel.csv",
                                           package = "MSstatsBioNet"))
    subnetwork <- suppressWarnings(suppressMessages(
        getSubnetworkFromIndra(input)))
    expect_named(subnetwork, c("nodes", "edges"))
})

test_that("save_network() and load_network() round-trip a network", {
    .mock_indra_with_clock()
    network <- .get_indra_network()
    file <- tempfile(fileext = ".rds")
    expect_identical(save_network(network, file), file)
    expect_identical(load_network(file), network)
})

test_that("save_network() warns when the network has no provenance", {
    .mock_indra_with_clock()
    network <- .get_indra_network()
    rebuilt <- list(nodes = network$nodes, edges = network$edges)
    file <- tempfile(fileext = ".rds")
    expect_warning(save_network(rebuilt, file), "no provenance")
    expect_true(file.exists(file))
    expect_no_warning(save_network(network, file))
})

test_that("save_network() checks its input", {
    .mock_indra_with_clock()
    network <- .get_indra_network()
    expect_error(save_network(network, "network.csv"), "ending in .rds")
    expect_error(save_network(network, c("a.rds", "b.rds")), "single path")
    bad <- network
    bad$edges$confidence <- 2
    file <- tempfile(fileext = ".rds")
    expect_error(save_network(bad, file), "confidence")
    expect_false(file.exists(file))
})

test_that("load_network() errors on a missing file or a network that fails the contract", {
    expect_error(load_network(file.path(tempdir(), "no-such-network.rds")),
                 "File not found")
    expect_error(load_network(NA_character_), "single path")
    file <- tempfile(fileext = ".rds")
    saveRDS(list(nodes = data.frame(id = "A"),
                 edges = data.frame(source = "A")), file)
    expect_error(load_network(file), "does not meet the current contract")
})
