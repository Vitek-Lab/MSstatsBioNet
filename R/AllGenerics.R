#' Get a network from a backend
#'
#' Sends the selected entities to a backend and returns the network that
#' answers the query, as \code{nodes} and \code{edges} tables.
#'
#' Only the rows of \code{entities} with \code{included_in_query = TRUE}
#' are sent to the backend. Every node the backend returns is matched
#' against all rows, so a node in the input gets \code{measured = TRUE} and
#' its statistics, and a node not in the input gets \code{measured = FALSE}
#' and \code{NA} statistics.
#'
#' \code{get_network()} prints the question it asks as a message, with the
#' number of entities, so the query in a saved script or log is readable
#' without the documentation.
#'
#' Confidence values are comparable within one backend, not across
#' backends. For INDRA, \code{confidence} is the INDRA belief score.
#'
#' @param backend a \code{NetworkBackend}, e.g. from
#' \code{\link{indra_backend}()}
#' @param entities entity table from \code{\link{prepare_entities}()},
#' grounded by \code{\link{convert_ids}()} and flagged by
#' \code{\link{select_entities}()}
#' @param query a \code{NetworkQuery} saying which question to ask, e.g.
#' \code{\link{subnetwork_query}()} (the default). See
#' \code{\link{network_queries}}.
#' @param interaction_types values of \code{edges$interaction} to keep
#' (INDRA statement types, e.g. \code{"Activation"}). \code{NULL} keeps all.
#' @param min_evidence minimum evidence count per edge
#' @param min_confidence minimum \code{confidence} per edge, in [0, 1].
#' Edges with no confidence score (\code{NA}) are dropped too, and a message
#' says how many. \code{NULL} applies no cutoff.
#' @param evidence_sources keeps edges with evidence from at least one of
#' these sources, e.g. \code{c("reach")}. \code{NULL} keeps all.
#' @param include_entities \code{"namespace:identifier"} groundings to add to
#' the query, e.g. \code{"HGNC:1234"}. Use this for entities outside the
#' input; to keep entities of the input that fail the cutoffs, use
#' \code{select_entities(force_include = )}.
#' @param ... passed to methods
#' @return list of \code{nodes} and \code{edges} data.frames that meets the
#' contract checked by \code{\link{validate_network}()}
#' @seealso \code{\link{network_queries}}, \code{\link{backend_capabilities}()}
#' @importFrom methods setGeneric setMethod
#' @export
#' @examples
#' input <- data.table::fread(system.file(
#'     "extdata/groupComparisonModel.csv",
#'     package = "MSstatsBioNet"
#' ))
#' indra <- indra_backend()
#' entities <- prepare_entities(input, entity_type = "protein",
#'                              id_type = "uniprot")
#' entities <- convert_ids(indra, entities)
#' entities <- select_entities(entities, pvalue_cutoff = 0.05)
#' network <- get_network(indra, entities, subnetwork_query(),
#'                        interaction_types = "Complex")
#' head(network$nodes)
#' head(network$edges)
setGeneric("get_network",
    function(backend, entities, query = subnetwork_query(),
             interaction_types = NULL, min_evidence = 1,
             min_confidence = NULL, evidence_sources = NULL,
             include_entities = NULL, ...)
        standardGeneric("get_network"),
    signature = c("backend", "query"))

# S4 dispatches a missing argument as class "missing", not on its default.
# The shared arguments are generic formals, so they don't reach `...` and
# must be passed on by name.
#' @rdname get_network
#' @export
setMethod("get_network", signature("NetworkBackend", "missing"),
    function(backend, entities, query, interaction_types = NULL,
             min_evidence = 1, min_confidence = NULL,
             evidence_sources = NULL, include_entities = NULL, ...) {
        get_network(backend, entities, subnetwork_query(),
                    interaction_types = interaction_types,
                    min_evidence = min_evidence,
                    min_confidence = min_confidence,
                    evidence_sources = evidence_sources,
                    include_entities = include_entities, ...)
    })

#' @rdname get_network
#' @export
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
#' left \code{NA}. \code{backend_capabilities(backend)$id_conversions} lists the
#' (entity type, identifier system) pairs a backend can convert.
#'
#' Ground the full results table, not only the significant rows, so that
#' \code{\link{get_network}()} can recognize every node that is in the
#' input.
#'
#' For INDRA, UniProt IDs and mnemonics are mapped through CoGEx, and gene
#' symbols and chemical names are grounded with Gilda. When an identifier
#' is a protein group (\code{"P1;P2"}) or a name with several candidates,
#' the groundings are \code{";"}-joined and positionally aligned.
#'
#' @param backend a \code{NetworkBackend}, e.g. from
#' \code{\link{indra_backend}()}
#' @param entities entity table from \code{\link{prepare_entities}()}
#' @param ... passed to methods
#' @return \code{entities} with the grounding columns filled in
#' @seealso \code{\link{get_entity_properties}()}
#' @export
#' @examples
#' df <- data.frame(Protein = c("P04637", "Q00610"))
#' entities <- prepare_entities(df, entity_type = "protein",
#'                              id_type = "uniprot")
#' convert_ids(indra_backend(), entities)
setGeneric("convert_ids",
    function(backend, entities, ...) standardGeneric("convert_ids"),
    signature = "backend")

#' @rdname convert_ids
#' @export
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
#' \code{backend_capabilities(backend)$entity_properties} lists the properties a
#' backend supports.
#'
#' For INDRA, the properties are \code{is_transcription_factor},
#' \code{is_kinase}, and \code{is_phosphatase}, looked up by gene symbol.
#' Only rows with a single HGNC grounding get values. A PTM site gets the
#' properties of its parent protein.
#'
#' @param backend a \code{NetworkBackend}, e.g. from
#' \code{\link{indra_backend}()}
#' @param entities entity table, grounded by \code{\link{convert_ids}()}
#' @param properties the properties to add. \code{NULL} adds every
#' property the backend supports.
#' @param ... passed to methods
#' @return \code{entities} with one column per property
#' @export
#' @examples
#' df <- data.frame(Protein = c("P04637", "Q00610"))
#' indra <- indra_backend()
#' entities <- prepare_entities(df, entity_type = "protein",
#'                              id_type = "uniprot")
#' entities <- convert_ids(indra, entities)
#' get_entity_properties(indra, entities, properties = "is_kinase")
setGeneric("get_entity_properties",
    function(backend, entities, properties = NULL, ...)
        standardGeneric("get_entity_properties"),
    signature = "backend")

#' @rdname get_entity_properties
#' @export
setMethod("get_entity_properties", "NetworkBackend",
    function(backend, entities, properties = NULL, ...) {
        stop(class(backend), " does not support get_entity_properties().",
             call. = FALSE)
    })

#' What a backend supports
#'
#' Lists the queries, identifier conversions, entity properties, and
#' limits of a backend, so a choice can be checked before a query is sent.
#'
#' @param backend a \code{NetworkBackend}, e.g. from
#' \code{\link{indra_backend}()}
#' @return named list:
#' \describe{
#'   \item{query_types}{the \code{query_type} values of the queries the
#'     backend answers, e.g. \code{"subnetwork"} for
#'     \code{\link{subnetwork_query}()}}
#'   \item{id_conversions}{for each \code{entity_type}, the \code{id_type}
#'     values \code{\link{convert_ids}()} can convert}
#'   \item{entity_properties}{the properties
#'     \code{\link{get_entity_properties}()} can add}
#'   \item{interaction_types}{the \code{edges$interaction} values the
#'     backend returns}
#'   \item{max_nodes}{for each query type, the largest number of
#'     groundings one query can take}
#' }
#' @seealso \code{\link{network_queries}}
#' @export
#' @examples
#' backend_capabilities(indra_backend())
setGeneric("backend_capabilities",
    function(backend) standardGeneric("backend_capabilities"),
    signature = "backend")

#' @rdname backend_capabilities
#' @export
setMethod("backend_capabilities", "NetworkBackend",
    function(backend) {
        stop(class(backend), " does not describe its backend_capabilities().",
             call. = FALSE)
    })
