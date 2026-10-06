# ----- get_evidence() and the backend_database registry (Phase 4b) -----

.make_evidence_edges <- function() {
    data.frame(
        source           = c("A", "B", "A"),
        target           = c("B", "C", "B"),
        interaction      = c("Activation", "Inhibition", "Activation"),
        site             = c(NA, "S473", "T308"),
        evidence_url     = c("https://example.com/1", "https://example.com/2",
                             "https://example.com/1"),
        statement_id     = c("111", "222", "111"),
        backend_database = "INDRA",
        stringsAsFactors = FALSE
    )
}

.make_indra_evidence <- function() {
    list(
        "111" = list(list(text = "A activates B.", pmid = "1001"),
                     list(text = "", pmid = "1002"),
                     list(text = "A binds B.")),
        "222" = list(list(text = "B inhibits C.", pmid = "2001"))
    )
}

# A second backend for the registry tests: its evidence names the backend.
setClass("TestEvidenceBackend", contains = "NetworkBackend",
         representation(label = "character"))
setMethod("get_evidence", "TestEvidenceBackend", function(backend, edges, ...) {
    data.frame(source = edges$source, target = edges$target,
               interaction = edges$interaction, site = edges$site,
               evidence_url = edges$evidence_url,
               statement_id = edges$statement_id,
               text = backend@label, pmid = "1",
               stringsAsFactors = FALSE)
})

describe("get_evidence() for INDRA", {

    test_that("copies each evidence sentence onto every edge with its hash", {
        local_mocked_bindings(
            .query_indra_evidence = function(stmt_hashes, cogex_url, ...) {
                .make_indra_evidence()
            }
        )
        evidence <- suppressMessages(capture.output(
            result <- get_evidence(indra_backend(), .make_evidence_edges())
        ))
        expect_equal(names(result),
                     c("source", "target", "interaction", "site",
                       "evidence_url", "statement_id", "text", "pmid"))
        # hash 111: two sentences with text, on edge rows 1 and 3; hash 222:
        # one sentence. The sentence with empty text is dropped.
        expect_equal(nrow(result), 5)
        expect_equal(result$text,
                     c("A activates B.", "A activates B.", "A binds B.",
                       "A binds B.", "B inhibits C."))
        expect_equal(result$site, c(NA, "T308", NA, "T308", "S473"))
        expect_equal(result$pmid, c("1001", "1001", "", "", "2001"))
    })

    test_that("asks for each statement hash once, at the backend's cogex_url", {
        requested <- NULL
        local_mocked_bindings(
            .query_indra_evidence = function(stmt_hashes, cogex_url, ...) {
                requested <<- list(hashes = stmt_hashes, url = cogex_url)
                .make_indra_evidence()
            }
        )
        backend <- indra_backend(cogex_url = "https://cogex.example.org")
        capture.output(get_evidence(backend, .make_evidence_edges()))
        expect_equal(requested$hashes, c("111", "222"))
        expect_equal(requested$url, "https://cogex.example.org")
    })

    test_that(".query_indra_evidence() posts to the evidence endpoint under cogex_url", {
        posted_url <- NULL
        local_mocked_bindings(
            POST = function(url, ...) {
                posted_url <<- url
                structure(list(), class = "response")
            },
            status_code = function(x) 200,
            content = function(x, ...) list()
        )
        capture.output(.query_indra_evidence("111", "https://cogex.example.org"))
        expect_equal(posted_url,
                     "https://cogex.example.org/api/get_evidences_for_stmt_hashes")
    })

    test_that("returns an empty table with a warning when there is no evidence", {
        local_mocked_bindings(
            .query_indra_evidence = function(stmt_hashes, cogex_url, ...) list()
        )
        capture.output(expect_warning(
            result <- get_evidence(indra_backend(), .make_evidence_edges()),
            "No evidence text found"
        ))
        expect_equal(result, .build_empty_evidence_table())
    })

    test_that("errors when edges lack a required column", {
        edges <- .make_evidence_edges()
        edges$evidence_url <- NULL
        expect_error(get_evidence(indra_backend(), edges),
                     "Missing required columns: evidence_url")
    })
})

test_that("get_evidence() errors for a backend without a method", {
    where <- environment()
    setClass("NoEvidenceBackend", contains = "NetworkBackend", where = where)
    on.exit(removeClass("NoEvidenceBackend", where = where), add = TRUE)
    expect_error(
        get_evidence(new("NoEvidenceBackend"), .make_evidence_edges()),
        "NoEvidenceBackend does not support get_evidence()", fixed = TRUE
    )
})

describe(".fetch_evidence() backend resolution", {

    test_that("uses the INDRA backend for INDRA edges", {
        used <- NULL
        local_mocked_bindings(
            .query_indra_evidence = function(stmt_hashes, cogex_url, ...) {
                used <<- cogex_url
                .make_indra_evidence()
            }
        )
        capture.output(result <- .fetch_evidence(.make_evidence_edges()))
        expect_equal(used, "https://discovery.indra.bio")
        expect_equal(nrow(result), 5)
    })

    test_that("splits a merged network by backend_database", {
        edges <- .make_evidence_edges()
        edges$backend_database <- c("First", "Second", "First")
        local_mocked_bindings(BACKEND_CONSTRUCTORS = list(
            First  = function() new("TestEvidenceBackend", label = "first"),
            Second = function() new("TestEvidenceBackend", label = "second")
        ))
        result <- .fetch_evidence(edges)
        expect_equal(result$statement_id, c("111", "111", "222"))
        expect_equal(result$text, c("first", "first", "second"))
    })

    test_that("an explicit backend takes every edge", {
        edges <- .make_evidence_edges()
        edges$backend_database <- c("INDRA", "Unknown", NA)
        result <- .fetch_evidence(
            edges, new("TestEvidenceBackend", label = "given"))
        expect_equal(result$text, rep("given", 3))
    })

    test_that("still resolves the backend after the network is rebuilt and subset", {
        network <- list(nodes = data.frame(id = c("A", "B", "C")),
                        edges = .make_evidence_edges())
        rebuilt <- list(nodes = network$nodes,
                        edges = network$edges[network$edges$source == "B", ])
        local_mocked_bindings(
            .query_indra_evidence = function(stmt_hashes, cogex_url, ...) {
                .make_indra_evidence()
            }
        )
        capture.output(result <- .fetch_evidence(rebuilt$edges))
        expect_equal(result$text, "B inhibits C.")
    })

    test_that("errors naming a backend_database with no default backend", {
        edges <- .make_evidence_edges()
        edges$backend_database <- c("INDRA", "BioGRID", "OmniPath")
        expect_error(.fetch_evidence(edges),
                     'No default backend for backend_database "BioGRID", "OmniPath"',
                     fixed = TRUE)
    })

    test_that("errors when backend_database is missing", {
        edges <- .make_evidence_edges()
        edges$backend_database <- NULL
        expect_error(.fetch_evidence(edges), "no `backend_database` column")
        edges$backend_database <- c("INDRA", NA, "INDRA")
        expect_error(.fetch_evidence(edges), "missing for some edges")
    })

    test_that("rejects a backend that is not a NetworkBackend", {
        expect_error(.fetch_evidence(.make_evidence_edges(), "INDRA"),
                     "`backend` must be NULL or a NetworkBackend")
    })

    test_that("returns an empty table for edges with no rows", {
        edges <- .make_evidence_edges()[0, ]
        expect_warning(result <- .fetch_evidence(edges), "No evidence text found")
        expect_equal(result, .build_empty_evidence_table())
    })
})

describe("downstream functions pass backend to get_evidence()", {

    given <- new("TestEvidenceBackend", label = "given")

    test_that("filterSubnetworkByContext() uses the given backend", {
        seen <- NULL
        local_mocked_bindings(
            .fetch_evidence = function(edges, backend = NULL) {
                seen <<- backend
                .build_empty_evidence_table()
            }
        )
        suppressWarnings(capture.output(filterSubnetworkByContext(
            data.frame(id = c("A", "B", "C")), .make_evidence_edges(),
            query = "kinase", backend = given)))
        expect_identical(seen, given)
    })

    test_that("the topic functions pass backend to the corpus", {
        seen <- list()
        local_mocked_bindings(
            .fetch_evidence = function(edges, backend = NULL) {
                seen[[length(seen) + 1]] <<- backend
                stop("stop after the evidence call")
            }
        )
        subnetwork <- list(nodes = data.frame(id = c("A", "B", "C")),
                           edges = .make_evidence_edges())
        expect_error(decomposeSubnetworkByTopic(subnetwork, n_topics = 2,
                                                backend = given),
                     "stop after")
        expect_error(compareTopicModels(subnetwork, seeds = 1:2, n_topics = 2,
                                        backend = given),
                     "stop after")
        expect_error(bootstrapTopicModels(subnetwork, n_boot = 2, n_topics = 2,
                                          backend = given),
                     "stop after")
        expect_error(decomposeSubnetworkIntoHierarchicalTopics(
            subnetwork, n_topics = 2, backend = given), "stop after")
        expect_length(seen, 4)
        for (backend in seen) expect_identical(backend, given)
    })

    test_that("the downstream functions reject an invalid backend before any query", {
        local_mocked_bindings(
            .fetch_evidence = function(...) stop("evidence was queried")
        )
        subnetwork <- list(nodes = data.frame(id = c("A", "B", "C")),
                           edges = .make_evidence_edges())
        message <- "`backend` must be NULL or a NetworkBackend"
        expect_error(filterSubnetworkByContext(
            subnetwork$nodes, subnetwork$edges, query = "kinase",
            backend = "INDRA"), message)
        expect_error(decomposeSubnetworkByTopic(subnetwork, backend = "INDRA"),
                     message)
        expect_error(compareTopicModels(subnetwork, seeds = 1:2,
                                        backend = "INDRA"), message)
        expect_error(bootstrapTopicModels(subnetwork, n_boot = 2,
                                          backend = "INDRA"), message)
        expect_error(decomposeSubnetworkIntoHierarchicalTopics(
            subnetwork, backend = "INDRA"), message)
    })
})
