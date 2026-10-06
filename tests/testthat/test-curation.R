# ----- get_curations() and filter_by_curation() (Phase 4c) -----

.make_curation_network <- function() {
    list(
        nodes = data.frame(id = c("A", "B", "C", "D", "E"),
                           stringsAsFactors = FALSE),
        edges = data.frame(
            source           = c("A", "B", "A", "C"),
            target           = c("B", "C", "C", "D"),
            interaction      = c("Activation", "Inhibition", "Complex",
                                 "Activation"),
            statement_id     = c("111", "222", "333", "111"),
            evidence_count   = c(5L, 2L, 1L, 3L),
            backend_database = "INDRA",
            stringsAsFactors = FALSE
        ),
        regulators = data.frame(id = "A")
    )
}

# A second backend for the registry tests: every statement has `incorrect`
# incorrect evidence.
setClass("TestCurationBackend", contains = "NetworkBackend",
         representation(incorrect = "integer"))
setMethod("get_curations", "TestCurationBackend", function(backend, edges, ...) {
    data.frame(statement_id = unique(edges$statement_id),
               incorrect_count = backend@incorrect,
               stringsAsFactors = FALSE)
})

test_that("indra_backend() defaults to the public INDRA database for curations", {
    expect_equal(indra_backend()@curation_url, "https://db.indra.bio")
})

test_that("indra_backend() rejects an invalid curation_url", {
    expect_error(indra_backend(curation_url = character(0)), "curation_url")
    expect_error(indra_backend(curation_url = ""), "curation_url")
    expect_error(indra_backend(curation_url = NA_character_), "curation_url")
})

describe("get_curations() for INDRA", {

    test_that("asks for each statement hash once, at the backend's curation_url", {
        requested <- list()
        local_mocked_bindings(
            .get_incorrect_curation_count = function(statement_id, curation_url) {
                requested[[length(requested) + 1]] <<- c(statement_id, curation_url)
                c("111" = 2, "222" = 0, "333" = 1)[[statement_id]]
            }
        )
        backend <- indra_backend(curation_url = "https://curation.example.org")
        result <- get_curations(backend, .make_curation_network()$edges)
        expect_equal(requested, list(
            c("111", "https://curation.example.org"),
            c("222", "https://curation.example.org"),
            c("333", "https://curation.example.org")))
        expect_equal(result, data.frame(statement_id = c("111", "222", "333"),
                                        incorrect_count = c(2L, 0L, 1L),
                                        stringsAsFactors = FALSE))
    })

    test_that("returns an empty table without a request for zero edges", {
        local_mocked_bindings(
            .get_incorrect_curation_count = function(...) stop("no request expected")
        )
        edges <- .make_curation_network()$edges[0, ]
        result <- get_curations(indra_backend(), edges)
        expect_equal(nrow(result), 0)
        expect_type(result$statement_id, "character")
        expect_type(result$incorrect_count, "integer")
    })

    test_that("errors when edges have no statement_id", {
        edges <- .make_curation_network()$edges
        edges$statement_id <- NULL
        expect_error(get_curations(indra_backend(), edges),
                     "Missing required columns: statement_id")
    })
})

describe(".get_incorrect_curation_count()", {

    .mock_curation_response <- function(env, curations, status = 200) {
        requested_url <- NULL
        local_mocked_bindings(
            GET = function(url, ...) {
                requested_url <<- url
                structure(list(), class = "response")
            },
            status_code = function(x) status,
            content = function(x, ...) "json",
            fromJSON = function(...) curations,
            .env = env
        )
        function() requested_url
    }

    test_that("counts each evidence curated as anything but correct once", {
        curations <- data.frame(
            tag = c("correct", "wrong_relation", "grounding", "wrong_relation"),
            source_hash = c(1, 2, 3, 2)
        )
        get_url <- .mock_curation_response(environment(), curations)
        count <- .get_incorrect_curation_count("-123", "https://curation.example.org")
        expect_equal(count, 2)
        expect_equal(get_url(), "https://curation.example.org/curation/list/-123")
    })

    test_that("counts 0 when a statement has no curations", {
        .mock_curation_response(environment(), list())
        expect_equal(.get_incorrect_curation_count("111"), 0)
    })

    test_that("warns and counts 0 when the request fails", {
        .mock_curation_response(environment(), list(), status = 500)
        expect_warning(count <- .get_incorrect_curation_count("111"),
                       "status code 500")
        expect_equal(count, 0)
    })
})

test_that("get_curations() errors for a backend without a method", {
    where <- environment()
    setClass("NoCurationBackend", contains = "NetworkBackend", where = where)
    on.exit(removeClass("NoCurationBackend", where = where), add = TRUE)
    expect_error(
        get_curations(new("NoCurationBackend"), .make_curation_network()$edges),
        "NoCurationBackend does not support get_curations()", fixed = TRUE
    )
})

describe("filter_by_curation()", {

    .mock_indra_curations <- function(env, counts = c("111" = 3, "222" = 0,
                                                     "333" = 1)) {
        local_mocked_bindings(
            .get_incorrect_curation_count = function(statement_id, curation_url) {
                counts[[statement_id]]
            },
            .env = env
        )
    }

    test_that("subtracts incorrect evidence and drops edges below min_evidence", {
        .mock_indra_curations(environment())
        expect_message(
            result <- filter_by_curation(.make_curation_network(),
                                         min_evidence = 2),
            "Dropping 2 edge(s) with fewer than 2 evidence", fixed = TRUE
        )
        # 111: 5 - 3 = 2 and 3 - 3 = 0; 222: 2 - 0 = 2; 333: 1 - 1 = 0
        expect_equal(result$edges$statement_id, c("111", "222"))
        expect_equal(result$edges$evidence_count, c(2L, 2L))
        expect_type(result$edges$evidence_count, "integer")
    })

    test_that("drops the nodes left without edges and keeps the rest", {
        .mock_indra_curations(environment())
        result <- suppressMessages(
            filter_by_curation(.make_curation_network(), min_evidence = 2))
        # D lost its only edge; E had no edge to begin with
        expect_equal(result$nodes$id, c("A", "B", "C", "E"))
        expect_equal(result$regulators, data.frame(id = "A"))
    })

    test_that("drops edges left with no evidence even when min_evidence is 0", {
        .mock_indra_curations(environment())
        result <- suppressMessages(
            filter_by_curation(.make_curation_network(), min_evidence = 0))
        expect_true(all(result$edges$evidence_count >= 1))
        expect_equal(result$edges$statement_id, c("111", "222"))
    })

    test_that("returns the network unchanged and silent when nothing is curated", {
        .mock_indra_curations(environment(),
                              counts = c("111" = 0, "222" = 0, "333" = 0))
        network <- .make_curation_network()
        expect_silent(result <- filter_by_curation(network))
        expect_equal(result, network)
    })

    test_that("uses the backend named in each edge's backend_database", {
        network <- .make_curation_network()
        network$edges$backend_database <- c("First", "Second", "First", "Second")
        local_mocked_bindings(BACKEND_CONSTRUCTORS = list(
            First  = function() new("TestCurationBackend", incorrect = 0L),
            Second = function() new("TestCurationBackend", incorrect = 1L)
        ))
        result <- suppressMessages(filter_by_curation(network))
        expect_equal(result$edges$statement_id, c("111", "222", "333", "111"))
        expect_equal(result$edges$evidence_count, c(5L, 1L, 1L, 2L))
    })

    test_that("an explicit backend takes every edge", {
        network <- .make_curation_network()
        network$edges$backend_database <- c("INDRA", "Unknown", NA, "INDRA")
        used_url <- NULL
        local_mocked_bindings(
            .get_incorrect_curation_count = function(statement_id, curation_url) {
                used_url <<- curation_url
                0
            }
        )
        result <- filter_by_curation(
            network, backend = indra_backend(curation_url = "https://db.example.org"))
        expect_equal(used_url, "https://db.example.org")
        expect_equal(nrow(result$edges), 4)
    })

    test_that("counts 0 for statements the backend returns no count for", {
        local_mocked_bindings(get_curations = function(backend, edges, ...) {
            data.frame(statement_id = "222", incorrect_count = 1L,
                       stringsAsFactors = FALSE)
        })
        result <- suppressMessages(filter_by_curation(.make_curation_network()))
        expect_equal(result$edges$evidence_count, c(5L, 1L, 1L, 3L))
    })

    test_that("errors when a backend returns a table that can't match the edges", {
        tables <- list(
            "a data.frame with columns statement_id and incorrect_count" =
                data.frame(statement_id = "111"),
            # numeric hashes lose precision above 2^53 and match no edge
            "a character statement_id with no NA" =
                data.frame(statement_id = c(111, 222, 333),
                           incorrect_count = 1L),
            "a character statement_id with no NA" =
                data.frame(statement_id = NA_character_, incorrect_count = 1L),
            "an incorrect_count of non-negative whole numbers" =
                data.frame(statement_id = "111", incorrect_count = -1L),
            "an incorrect_count of non-negative whole numbers" =
                data.frame(statement_id = "111", incorrect_count = 0.5),
            "an incorrect_count of non-negative whole numbers" =
                data.frame(statement_id = "111", incorrect_count = NA_integer_),
            "an incorrect_count of non-negative whole numbers" =
                data.frame(statement_id = "111", incorrect_count = "1")
        )
        for (i in seq_along(tables)) {
            local_mocked_bindings(get_curations = function(backend, edges, ...) {
                tables[[i]]
            })
            expect_error(filter_by_curation(.make_curation_network()),
                         paste0("get_curations() for IndraBackend must return ",
                                names(tables)[i]), fixed = TRUE)
        }
    })

    test_that("errors for a backend_database with no default backend", {
        network <- .make_curation_network()
        network$edges$backend_database <- "OmniPath"
        expect_error(filter_by_curation(network),
                     'No default backend for backend_database "OmniPath"',
                     fixed = TRUE)
    })

    test_that("errors without backend_database unless a backend is given", {
        network <- .make_curation_network()
        network$edges$backend_database <- NULL
        expect_error(filter_by_curation(network), "no `backend_database` column")
    })

    test_that("checks its input before any request", {
        local_mocked_bindings(
            .get_incorrect_curation_count = function(...) stop("no request expected")
        )
        network <- .make_curation_network()
        expect_error(filter_by_curation(network$edges), "must be a list")
        no_count <- network
        no_count$edges$evidence_count <- NULL
        expect_error(filter_by_curation(no_count), "evidence_count")
        no_id <- network
        no_id$nodes$id <- NULL
        expect_error(filter_by_curation(no_id), "nodes` is missing the column: id")
        expect_error(filter_by_curation(network, min_evidence = "2"),
                     "single number")
        expect_error(filter_by_curation(network, min_evidence = c(1, 2)),
                     "single number")
        expect_error(filter_by_curation(network, min_evidence = NA),
                     "single number")
        expect_error(filter_by_curation(network, backend = "INDRA"),
                     "must be NULL or a NetworkBackend")
    })
})
