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
#' @param interaction_types values of \code{edges$interaction} to keep
#' (INDRA statement types, e.g. \code{"Activation"}). \code{NULL} keeps all.
#' @param min_evidence minimum evidence count per edge
#' @param evidence_sources keeps edges with evidence from at least one of
#' these sources, e.g. \code{c("reach")}. \code{NULL} keeps all.
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
             interaction_types = NULL, min_evidence = 1,
             evidence_sources = NULL, include_entities = NULL, ...)
        standardGeneric("get_network"),
    signature = c("backend", "query"))

# S4 dispatches a missing argument as class "missing", not on its default.
# The shared arguments are generic formals, so they don't reach `...` and
# must be passed on by name.
setMethod("get_network", signature("NetworkBackend", "missing"),
    function(backend, entities, query, interaction_types = NULL,
             min_evidence = 1, evidence_sources = NULL,
             include_entities = NULL, ...) {
        get_network(backend, entities, subnetwork_query(),
                    interaction_types = interaction_types,
                    min_evidence = min_evidence,
                    evidence_sources = evidence_sources,
                    include_entities = include_entities, ...)
    })

setMethod("get_network", signature("NetworkBackend", "NetworkQuery"),
    function(backend, entities, query, ...) {
        stop(class(backend), " does not support ", class(query), ".",
             call. = FALSE)
    })

#' Ground entities in a backend's namespaces
#'
#' Fills in \code{namespace}, \code{entity_id}, and \code{entity_name} from
#' each row's \code{id} (\code{parent_id} for \code{ptm_site} rows), read
#' as the identifier system in \code{id_type}. Rows that don't ground are
#' left \code{NA}.
#'
#' Internal until the end of Phase 3 of the API refactor.
#'
#' @param backend a \code{NetworkBackend}, e.g. from \code{indra_backend()}
#' @param entities entity table from \code{prepare_entities()}
#' @param ... passed to methods
#' @return \code{entities} with the grounding columns filled in
#' @keywords internal
#' @noRd
setGeneric("convert_ids",
    function(backend, entities, ...) standardGeneric("convert_ids"),
    signature = "backend")

setMethod("convert_ids", "NetworkBackend",
    function(backend, entities, ...) {
        stop(class(backend), " does not support convert_ids().",
             call. = FALSE)
    })

#' Add a backend's properties of each entity
#'
#' Adds one column per property, e.g. \code{is_kinase}. A property is
#' \code{NA} for rows whose \code{entity_type} it doesn't apply to, and for
#' rows the backend has no answer for.
#'
#' Internal until the end of Phase 3 of the API refactor.
#'
#' @param backend a \code{NetworkBackend}, e.g. from \code{indra_backend()}
#' @param entities entity table, grounded by \code{convert_ids()}
#' @param properties the properties to add. \code{NULL} adds every
#' property the backend supports.
#' @param ... passed to methods
#' @return \code{entities} with one column per property
#' @keywords internal
#' @noRd
setGeneric("get_entity_properties",
    function(backend, entities, properties = NULL, ...)
        standardGeneric("get_entity_properties"),
    signature = "backend")

setMethod("get_entity_properties", "NetworkBackend",
    function(backend, entities, properties = NULL, ...) {
        stop(class(backend), " does not support get_entity_properties().",
             call. = FALSE)
    })
