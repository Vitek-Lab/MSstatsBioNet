# ----- INDRA backend and get_network() (Phase 2 of the API refactor) -----

.selected_input <- function() {
    input <- data.table::fread(
        system.file("extdata/groupComparisonModel.csv", package = "MSstatsBioNet")
    )
    .filterGetSubnetworkFromIndraInput(input, pvalueCutoff = NULL,
                                       logfc_cutoff = NULL,
                                       force_include_other = NULL,
                                       include_infinite_fc = FALSE,
                                       direction = "both")
}

.mock_indra_response <- function(env = parent.frame()) {
    local_mocked_bindings(
        .callIndraCogexApi = function(ns, ids, fio, cogex_url) {
            readRDS(system.file("extdata/indraResponse.rds", package = "MSstatsBioNet"))
        },
        .env = env
    )
}

test_that("indra_backend() defaults to the public CoGEx URL", {
    backend <- indra_backend()
    expect_s4_class(backend, "IndraBackend")
    expect_s4_class(backend, "NetworkBackend")
    expect_equal(backend@cogex_url, "https://discovery.indra.bio")
})

test_that("indra_backend() rejects an invalid cogex_url", {
    expect_error(indra_backend(cogex_url = character(0)), "cogex_url")
    expect_error(indra_backend(cogex_url = c("a", "b")), "cogex_url")
    expect_error(indra_backend(cogex_url = ""), "cogex_url")
    expect_error(indra_backend(cogex_url = NA_character_), "cogex_url")
})

test_that("get_network() reproduces the pinned golden output", {
    .mock_indra_response()
    network <- get_network(indra_backend(), .selected_input(), subnetwork_query())
    golden <- readRDS(test_path("_fixtures", "golden_subnetwork.rds"))
    expect_identical(network, golden)
})

test_that("get_network() runs subnetwork_query() when query is missing", {
    .mock_indra_response()
    input <- .selected_input()
    expect_identical(
        get_network(indra_backend(), input),
        get_network(indra_backend(), input, subnetwork_query())
    )
})

test_that("get_network() passes the arguments on to the old filters", {
    .mock_indra_response()
    input <- .selected_input()
    network <- get_network(indra_backend(), input,
                           statement_types = "Complex", min_evidence = 2)
    expect_true(all(network$edges$interaction == "Complex"))
    expect_true(all(network$edges$evidence_count >= 2))

    network <- get_network(indra_backend(), input,
                           statement_types = c("Activation", "Phosphorylation"))
    expect_equal(nrow(network$nodes), 3)
    expect_equal(nrow(network$edges), 2)
})

test_that("get_network() sends the backend's cogex_url", {
    sent_url <- NULL
    local_mocked_bindings(
        .callIndraCogexApi = function(ns, ids, fio, cogex_url) {
            sent_url <<- cogex_url
            readRDS(system.file("extdata/indraResponse.rds", package = "MSstatsBioNet"))
        }
    )
    get_network(indra_backend(cogex_url = "https://cogex.example.org"),
                .selected_input())
    expect_equal(sent_url, "https://cogex.example.org")
})

test_that(".callIndraCogexApi() posts to the subnetwork endpoint under cogex_url", {
    posted_url <- NULL
    local_mocked_bindings(
        POST = function(url, ...) {
            posted_url <<- url
            structure(list(), class = "response")
        },
        content = function(x, ...) list()
    )
    .callIndraCogexApi("HGNC", "1925", NULL, "https://cogex.example.org")
    expect_equal(posted_url,
                 "https://cogex.example.org/api/indra_subnetwork_relations")
})

test_that("get_network() errors for an unsupported backend and query pair", {
    where <- environment()
    setClass("UnsupportedQuery", contains = "NetworkQuery", where = where)
    on.exit(removeClass("UnsupportedQuery", where = where), add = TRUE)
    expect_error(
        get_network(indra_backend(), .selected_input(),
                    new("UnsupportedQuery")),
        "IndraBackend does not support UnsupportedQuery"
    )
})

test_that("get_network() validates its input before calling INDRA", {
    local_mocked_bindings(
        .callIndraCogexApi = function(ns, ids, fio, cogex_url) {
            stop("INDRA should not be called")
        }
    )
    input <- .selected_input()
    expect_error(get_network(indra_backend(), input[0, ]),
                 "at least one protein")
    no_entity_id <- as.data.frame(input)
    no_entity_id$EntityId <- NULL
    expect_error(get_network(indra_backend(), no_entity_id),
                 "missing required column\\(s\\): EntityId")
    expect_error(get_network(indra_backend(), input, sources = 1),
                 "sources_filter must be a character vector")
})
