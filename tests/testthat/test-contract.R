.valid_network <- function() {
    list(
        nodes = data.frame(
            id = c("A", "B", "C"),
            entity_name = c("CHEK1", "CDC25A", "TP53"),
            log2FC = c(1.2, -0.5, NA),
            adj.pvalue = c(0.01, 0.2, NA),
            stringsAsFactors = FALSE
        ),
        edges = data.frame(
            source = c("A", "B"),
            target = c("B", "C"),
            interaction = c("Phosphorylation", "Complex"),
            directed = c(TRUE, FALSE),
            site = c("S76", NA),
            confidence = c(0.99, NA),
            evidence_count = c(12L, 1L),
            provenance_url = c(
                "https://db.indra.bio/statements/from_hash/1?format=html",
                "https://db.indra.bio/statements/from_hash/-2?format=html"),
            statement_id = c("1", "-2"),
            source_db = "INDRA",
            query_type = "subnetwork",
            stringsAsFactors = FALSE
        )
    )
}

test_that("validate_network accepts a network that meets the contract", {
    network <- .valid_network()
    expect_identical(validate_network(network), network)
    expect_invisible(validate_network(network))
})

test_that("validate_network accepts empty edges with the right columns", {
    network <- .valid_network()
    network$edges <- network$edges[0, ]
    expect_silent(validate_network(network))
})

test_that("validate_network requires a list with nodes and edges data.frames", {
    expect_error(validate_network(data.frame()), "list with 'nodes' and 'edges'")
    expect_error(validate_network(list(nodes = data.frame(id = "A"), edges = 1)),
                 "must be data.frames")
})

test_that("validate_network reports missing required edge columns", {
    network <- .valid_network()
    network$edges$confidence <- NULL
    network$edges$source_db <- NULL
    expect_error(validate_network(network),
                 "edges is missing required column\\(s\\): confidence, source_db")
})

test_that("validate_network reports a missing nodes$id column", {
    network <- .valid_network()
    network$nodes$id <- NULL
    expect_error(validate_network(network), "nodes is missing required column")
})

test_that("validate_network checks column types", {
    network <- .valid_network()
    network$edges$statement_id <- c(1, -2)
    expect_error(validate_network(network), "edges\\$statement_id must be character")

    network <- .valid_network()
    network$edges$directed <- c("yes", "no")
    expect_error(validate_network(network), "edges\\$directed must be logical")

    network <- .valid_network()
    network$edges$evidence_count <- c(1.5, 2)
    expect_error(validate_network(network), "edges\\$evidence_count must be integer")

    network <- .valid_network()
    network$nodes$log2FC <- c("1", "2", NA)
    expect_error(validate_network(network), "nodes\\$log2FC must be numeric")
})

test_that("validate_network accepts whole-number doubles for evidence_count", {
    network <- .valid_network()
    network$edges$evidence_count <- c(12, 1)
    expect_silent(validate_network(network))
})

test_that("validate_network accepts all-NA logical columns of any type", {
    network <- .valid_network()
    network$edges$site <- NA
    network$edges$confidence <- NA
    expect_silent(validate_network(network))
})

test_that("validate_network checks the statement-type vocabulary", {
    network <- .valid_network()
    network$edges$interaction[1] <- "Phosphorylated"
    expect_error(validate_network(network),
                 "outside the statement-type vocabulary: Phosphorylated")
})

test_that("validate_network rejects NA in identifying edge columns", {
    network <- .valid_network()
    network$edges$statement_id[1] <- NA
    expect_error(validate_network(network), "edges\\$statement_id must not be NA")
})

test_that("validate_network checks confidence and evidence_count ranges", {
    network <- .valid_network()
    network$edges$confidence[1] <- 1.5
    expect_error(validate_network(network), "confidence must be in \\[0, 1\\]")

    network <- .valid_network()
    network$edges$evidence_count[1] <- 0L
    expect_error(validate_network(network), "evidence_count must be at least 1")
})

test_that("validate_network rejects NA or infinite evidence_count", {
    network <- .valid_network()
    network$edges$evidence_count <- c(12, NA)
    expect_error(validate_network(network),
                 "evidence_count must not be NA or infinite")

    network$edges$evidence_count <- c(12, Inf)
    expect_error(validate_network(network),
                 "evidence_count must not be NA or infinite")

    # An all-NA logical column passes the type check, but not the NA check
    network$edges$evidence_count <- NA
    expect_error(validate_network(network),
                 "evidence_count must not be NA or infinite")

    # A character column gets both the type error and the NA error
    network$edges$evidence_count <- c("12", NA)
    err <- tryCatch(validate_network(network),
                    error = function(e) conditionMessage(e))
    expect_match(err, "edges\\$evidence_count must be integer")
    expect_match(err, "evidence_count must not be NA or infinite")
})

test_that("validate_network requires directed == FALSE for symmetric types", {
    network <- .valid_network()
    network$edges$directed[2] <- TRUE   # the Complex edge
    expect_error(validate_network(network),
                 "directed must be FALSE for symmetric statement types")

    network <- .valid_network()
    network$edges$interaction[2] <- "Association"
    network$edges$directed[2] <- TRUE
    expect_error(validate_network(network), "1 row\\(s\\) have directed == TRUE")

    # Other statement types may be undirected
    network <- .valid_network()
    network$edges$directed[1] <- FALSE   # the Phosphorylation edge
    expect_silent(validate_network(network))

    # NA directed is reported by the NA check only
    network <- .valid_network()
    network$edges$directed[2] <- NA
    err <- tryCatch(validate_network(network),
                    error = function(e) conditionMessage(e))
    expect_match(err, "edges\\$directed must not be NA")
    expect_no_match(err, "symmetric statement types")
})

test_that("validate_network checks the site format", {
    network <- .valid_network()
    network$edges$site[1] <- "S76;T80"
    expect_silent(validate_network(network))

    network$edges$site[1] <- "Ser76"
    expect_error(validate_network(network), "edges\\$site must look like 'S148'")
})

test_that("validate_network requires every edge endpoint to be a node", {
    network <- .valid_network()
    network$edges$target[2] <- "D"
    expect_error(validate_network(network), "not found in nodes\\$id: D")
})

test_that("validate_network rejects NA node ids", {
    network <- .valid_network()
    network$nodes$id[3] <- NA
    network$edges <- network$edges[1, ]
    expect_error(validate_network(network), "nodes\\$id must not be NA")
})

test_that("validate_network accepts one node row per PTM site of a protein", {
    network <- .valid_network()
    network$nodes <- rbind(network$nodes, network$nodes[1, ])
    network$nodes$site <- c("S317", NA, NA, "S345")
    expect_silent(validate_network(network))
})

test_that("validate_network checks node vocabularies when the columns exist", {
    network <- .valid_network()
    network$nodes$entity_type <- c("protein", "protein", "enzyme")
    expect_error(validate_network(network),
                 "nodes\\$entity_type has value\\(s\\) outside the vocabulary: enzyme")

    network <- .valid_network()
    network$nodes$node_role <- c("query", "query", "regulator")
    expect_silent(validate_network(network))
})

test_that("validate_network requires NA statistics on latent nodes", {
    network <- .valid_network()
    network$nodes$measured <- c(TRUE, TRUE, FALSE)
    expect_silent(validate_network(network))

    network$nodes$log2FC[3] <- 0
    network$nodes$adj.pvalue[3] <- 1
    expect_error(validate_network(network),
                 "nodes\\$log2FC must be NA for nodes with measured == FALSE")
    expect_error(validate_network(network),
                 "nodes\\$adj.pvalue must be NA for nodes with measured == FALSE")
})

test_that("validate_network reports a non-logical measured column", {
    network <- .valid_network()
    network$nodes$measured <- c("yes", "yes", "no")
    network$nodes$log2FC[3] <- 0
    err <- tryCatch(validate_network(network),
                    error = function(e) conditionMessage(e))
    expect_match(err, "nodes\\$measured must be logical")
    expect_no_match(err, "must be NA for nodes with measured == FALSE")
})

test_that("validate_network lists every problem in one error", {
    network <- .valid_network()
    network$edges$interaction[1] <- "Unknown"
    network$edges$confidence[2] <- -1
    err <- tryCatch(validate_network(network), error = function(e) conditionMessage(e))
    expect_match(err, "statement-type vocabulary")
    expect_match(err, "confidence must be in")
})

test_that("UNDIRECTED_STATEMENT_TYPES are part of the statement vocabulary", {
    expect_true(all(MSstatsBioNet:::UNDIRECTED_STATEMENT_TYPES %in%
                        MSstatsBioNet:::STATEMENT_TYPES))
})
