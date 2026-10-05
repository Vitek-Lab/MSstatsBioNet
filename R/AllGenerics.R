#' Get a network from a backend
#'
#' Dispatches on \code{backend} and \code{query}. A missing \code{query}
#' runs \code{subnetwork_query()}.
#'
#' Internal until the entity model is added (Phase 3 of the API refactor).
#' Until then \code{entities} is the groupComparison table annotated by
#' \code{annotateProteinInfoFromIndra()} and already filtered to the
#' selected rows.
#'
#' @param backend a \code{NetworkBackend}, e.g. from \code{indra_backend()}
#' @param entities annotated groupComparison table of the selected rows
#' @param query a \code{NetworkQuery}, e.g. from \code{subnetwork_query()}
#' @param statement_types statement types to keep. \code{NULL} keeps all.
#' @param min_evidence minimum evidence count per edge
#' @param sources evidence sources to keep, e.g. \code{c("reach")}.
#' \code{NULL} keeps all.
#' @param include_entities \code{"namespace:identifier"} strings to add to
#' the query, e.g. \code{"HGNC:1234"}
#' @param ... passed to methods
#' @return list of \code{nodes} and \code{edges} that meets the contract
#' checked by \code{validate_network()}
#' @importFrom methods setGeneric setMethod
#' @keywords internal
#' @noRd
setGeneric("get_network",
    function(backend, entities, query = subnetwork_query(),
             statement_types = NULL, min_evidence = 1, sources = NULL,
             include_entities = NULL, ...)
        standardGeneric("get_network"),
    signature = c("backend", "query"))

# S4 dispatches a missing argument as class "missing", not on its default.
# The shared arguments are generic formals, so they don't reach `...` and
# must be passed on by name.
setMethod("get_network", signature("NetworkBackend", "missing"),
    function(backend, entities, query, statement_types = NULL,
             min_evidence = 1, sources = NULL, include_entities = NULL,
             ...) {
        get_network(backend, entities, subnetwork_query(),
                    statement_types = statement_types,
                    min_evidence = min_evidence, sources = sources,
                    include_entities = include_entities, ...)
    })

setMethod("get_network", signature("NetworkBackend", "NetworkQuery"),
    function(backend, entities, query, ...) {
        stop(class(backend), " does not support ", class(query), ".",
             call. = FALSE)
    })
