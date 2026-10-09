# Shared contract checks for every network backend. Each backend's tests
# call expect_backend_contract() with its HTTP layer mocked.

# Query constructor for each query_type a backend can list in
# backend_capabilities()$query_types
.contract_query_constructors <- list(
    subnetwork = subnetwork_query
)

# Checks backend_capabilities(backend), then runs get_network() for every
# query type it lists and checks the network: it meets the contract, has
# edges from one backend_database, from the query type asked, and of the
# interaction types the backend lists, its backend_database has a default
# backend for get_evidence() and get_curations(), and merging it with
# itself gives it back.
expect_backend_contract <- function(backend, entities, ...) {
    capabilities <- backend_capabilities(backend)
    expect_type(capabilities, "list")
    expect_true(all(c("query_types", "id_conversions", "interaction_types",
                      "max_nodes") %in% names(capabilities)))
    expect_type(capabilities$query_types, "character")
    expect_true(all(capabilities$interaction_types %in% INTERACTION_TYPES))
    expect_true(all(names(capabilities$id_conversions) %in% ENTITY_TYPES))
    if (!is.null(capabilities$evidence_sources)) {
        expect_type(capabilities$evidence_sources, "list")
        expect_true(all(c("database", "text_mined") %in%
                        names(capabilities$evidence_sources)))
        expect_true(all(vapply(capabilities$evidence_sources, is.character,
                               logical(1))))
    }

    for (query_type in capabilities$query_types) {
        constructor <- .contract_query_constructors[[query_type]]
        expect_false(is.null(constructor),
                     info = paste("no query constructor for", query_type))
        network <- suppressMessages(
            get_network(backend, entities, constructor(), ...))
        expect_no_error(validate_network(network))
        expect_gt(nrow(network$edges), 0)
        database <- unique(network$edges$backend_database)
        expect_length(database, 1)
        expect_true(database %in% names(BACKEND_CONSTRUCTORS))
        expect_true(all(network$edges$query_type == query_type))
        expect_true(all(network$edges$interaction %in%
                        capabilities$interaction_types))
        merged <- merge_networks(network, network)
        expect_identical(merged$edges, .reset_rownames(network$edges))
        expect_identical(merged$nodes, .reset_rownames(network$nodes))
    }
    invisible(backend)
}

.reset_rownames <- function(table) {
    rownames(table) <- NULL
    table
}
