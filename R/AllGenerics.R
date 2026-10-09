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
#' Each node is one row of \code{entities} (a PTM site row is drawn on its
#' parent protein's node). When several rows share a grounding, e.g. two
#' isoforms that both ground to the same gene, or a protein and a protein
#' group that contains it, each of their nodes gets the backend's edges,
#' with its own statistics. So one backend statement can give several
#' edges, and edges can share a \code{statement_id}. A statement from a
#' grounding to itself, such as a homodimer, gives each matching node a
#' self-loop, and no edges between those nodes. To count statements rather
#' than edges, count unique \code{backend_database} and
#' \code{statement_id} pairs.
#'
#' \code{get_network()} prints the question it asks as a message, with the
#' number of entities, so the query in a saved script or log is readable
#' without the documentation.
#'
#' Confidence values are comparable within one backend, not across
#' backends. For INDRA, \code{confidence} is the INDRA belief score; for
#' STRING, the combined score.
#'
#' Each backend reads only its own groundings of \code{entities}, so one
#' table grounded by several backends' \code{\link{convert_ids}()} serves
#' all of them, and their networks share node IDs for
#' \code{\link{merge_networks}()}.
#'
#' @param backend a \code{NetworkBackend}, e.g. from
#' \code{\link{indra_backend}()} or \code{\link{string_backend}()}
#' @param entities entity table from \code{\link{prepare_entities}()},
#' grounded by \code{\link{convert_ids}()} and flagged by
#' \code{\link{select_entities}()}
#' @param query a \code{NetworkQuery} saying which question to ask, e.g.
#' \code{\link{subnetwork_query}()} (the default). See
#' \code{\link{network_queries}}.
#' @param interaction_types values of \code{edges$interaction} to keep
#' (INDRA statement types, e.g. \code{"Activation"}). \code{NULL} keeps all.
#' @param min_evidence minimum evidence count per edge. Edges with no
#' evidence count (\code{NA}, as from STRING) pass the default of 1, and
#' are dropped by a higher cutoff, with a message saying how many.
#' @param min_confidence minimum \code{confidence} per edge, in [0, 1].
#' Edges with no confidence score (\code{NA}) are dropped too, and a message
#' says how many. \code{NULL} applies no cutoff. STRING gets the cutoff as
#' its \code{required_score}; STRING's website uses 0.4 by default.
#' @param evidence_sources keeps edges with evidence from at least one of
#' these sources, e.g. \code{c("reach")}. \code{NULL} keeps all.
#' \code{backend_capabilities(backend)$evidence_sources} lists a backend's
#' sources; for INDRA, \code{\link{INDRA_DATABASE_SOURCES}} keeps edges
#' with curated-database evidence and
#' \code{\link{INDRA_TEXT_MINED_SOURCES}} edges with text-mined evidence.
#' For STRING, the sources are its evidence channels, e.g.
#' \code{"experiments"}, and an edge is kept when one of them has a
#' non-zero score.
#' @param include_entities \code{"namespace:identifier"} groundings to add to
#' the query, e.g. \code{"HGNC:1234"}, or for STRING
#' \code{"STRING:9606.ENSP00000269305"}. Use this for entities outside the
#' input; to keep entities of the input that fail the cutoffs, use
#' \code{select_entities(force_include = )}.
#' @param ... passed to methods
#' @return list of \code{nodes} and \code{edges} data.frames that meets the
#' contract checked by \code{\link{validate_network}()}, and
#' \code{provenance}: a data.frame with one row recording the query, with
#' the columns \code{backend_database}, \code{query_type},
#' \code{retrieved_at} (when the response arrived, in UTC),
#' \code{backend_version} (e.g. \code{"12.5"} for STRING; \code{NA} for
#' INDRA, which has no data versions), \code{backend_url}, \code{organism} (NCBI taxon IDs of the
#' queried rows), \code{parameters} (the query arguments, as JSON), and
#' \code{package_version} (of MSstatsBioNet). It joins to the edges on
#' \code{backend_database} and \code{query_type}.
#' \code{\link{merge_networks}()} combines the provenance of its networks,
#' and \code{\link{save_network}()} saves it with the network. It is lost
#' when the network is rebuilt with \code{list(nodes = , edges = )}.
#' @seealso \code{\link{network_queries}}, \code{\link{backend_capabilities}()}
#' @importFrom methods setGeneric setMethod
#' @export
#' @examples
#' input <- data.table::fread(system.file(
#'     "extdata/groupComparisonModel.csv",
#'     package = "MSstatsBioNet"
#' ))
#' entities <- prepare_entities(input, entity_type = "protein",
#'                              id_type = "uniprot")
#' \donttest{
#' indra <- indra_backend()
#' entities <- convert_ids(indra, entities)
#' entities <- select_entities(entities, pvalue_cutoff = 0.05)
#' network <- get_network(indra, entities, subnetwork_query(),
#'                        interaction_types = "Complex")
#' head(network$nodes)
#' head(network$edges)
#' }
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
#' A backend replaces only the groundings in its own namespaces and keeps
#' the others, so a table grounded by \code{indra_backend()} and then by
#' \code{string_backend()} has both groundings, e.g. \code{namespace =
#' "HGNC;STRING"}, and serves both backends.
#'
#' For STRING, only UniProt accessions are converted (\code{id_type =
#' "uniprot"}), each to one STRING protein of the row's organism, and the
#' namespace is \code{"STRING"}. An isoform (\code{"P04637-2"}) gets the
#' protein of its canonical accession. Identifiers that aren't accessions
#' give an error before any request, since STRING's search maps other text
#' to the wrong protein, and accessions STRING doesn't know are listed in a
#' message.
#'
#' For INDRA, UniProt IDs and mnemonics are mapped through CoGEx, and gene
#' symbols and chemical names are grounded with Gilda. When an identifier
#' is a protein group (\code{"P1;P2"}) or a name with several candidates,
#' the groundings are \code{";"}-joined and positionally aligned.
#'
#' @param backend a \code{NetworkBackend}, e.g. from
#' \code{\link{indra_backend}()} or \code{\link{string_backend}()}
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

#' Get the evidence behind network edges
#'
#' Looks up the evidence sentences and their PubMed IDs for each edge,
#' using the edge's \code{statement_id}. Edges that share a
#' \code{statement_id} get the same evidence.
#'
#' For INDRA, the evidence comes from CoGEx. Evidence without text is left
#' out, and \code{pmid} is \code{""} for evidence from a source with no
#' PubMed ID.
#'
#' \code{\link{filterSubnetworkByContext}()} and the topic functions call
#' \code{get_evidence()} with the backend named in each edge's
#' \code{backend_database}, unless their \code{backend} argument is given.
#'
#' @param backend a \code{NetworkBackend}, e.g. from
#' \code{\link{indra_backend}()}
#' @param edges the \code{edges} of a network from
#' \code{\link{get_network}()}. Needs the columns \code{source},
#' \code{target}, \code{interaction}, \code{site}, \code{evidence_url}, and
#' \code{statement_id}.
#' @param ... passed to methods
#' @return data.frame with one row per (edge, evidence sentence) pair and
#' the columns \code{source}, \code{target}, \code{interaction},
#' \code{site}, \code{evidence_url}, \code{statement_id}, \code{text}, and
#' \code{pmid}. It has no rows, with a warning, when no edge has evidence
#' text.
#' @seealso \code{\link{filterSubnetworkByContext}()}
#' @export
#' @examples
#' \donttest{
#' input <- data.table::fread(system.file(
#'     "extdata/groupComparisonModel.csv",
#'     package = "MSstatsBioNet"
#' ))
#' indra <- indra_backend()
#' entities <- prepare_entities(input, entity_type = "protein",
#'                              id_type = "uniprot")
#' entities <- convert_ids(indra, entities)
#' entities <- select_entities(entities, pvalue_cutoff = 0.05)
#' network <- get_network(indra, entities, interaction_types = "Complex")
#' evidence <- get_evidence(indra, head(network$edges, 2))
#' head(evidence[, c("source", "target", "pmid", "text")])
#' }
setGeneric("get_evidence",
    function(backend, edges, ...) standardGeneric("get_evidence"),
    signature = "backend")

#' @rdname get_evidence
#' @export
setMethod("get_evidence", "NetworkBackend",
    function(backend, edges, ...) {
        stop(class(backend), " does not support get_evidence().",
             call. = FALSE)
    })

#' Get the curations of edges from a backend
#'
#' Curators mark single pieces of evidence for a statement as correct or
#' incorrect. \code{get_curations()} counts, for each statement, the
#' evidence curated as incorrect. INDRA curations come from the INDRA
#' database, with one request per statement. Most users call
#' \code{\link{filter_by_curation}()}, which subtracts these counts from
#' \code{evidence_count}.
#'
#' @inheritParams get_evidence
#' @param edges the \code{edges} of a network from
#' \code{\link{get_network}()}. Needs the column \code{statement_id}.
#' @return data.frame with one row per unique \code{statement_id} and the
#' columns \code{statement_id} (character) and \code{incorrect_count}
#' (integer). A failed request warns and counts as 0. A method for another
#' backend may leave out statements with no curations, which
#' \code{\link{filter_by_curation}()} counts as 0, but must return
#' \code{statement_id} as character: \code{filter_by_curation()} errors
#' otherwise, since a numeric hash can lose precision and match no edge.
#' @seealso \code{\link{filter_by_curation}()}
#' @export
#' @examples
#' \donttest{
#' input <- data.table::fread(system.file(
#'     "extdata/groupComparisonModel.csv",
#'     package = "MSstatsBioNet"
#' ))
#' indra <- indra_backend()
#' entities <- prepare_entities(input, entity_type = "protein",
#'                              id_type = "uniprot")
#' entities <- convert_ids(indra, entities)
#' entities <- select_entities(entities, pvalue_cutoff = 0.05)
#' network <- get_network(indra, entities, interaction_types = "Complex")
#' get_curations(indra, head(network$edges, 2))
#' }
setGeneric("get_curations",
    function(backend, edges, ...) standardGeneric("get_curations"),
    signature = "backend")

#' @rdname get_curations
#' @export
setMethod("get_curations", "NetworkBackend",
    function(backend, edges, ...) {
        stop(class(backend), " does not support get_curations().",
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
#'   \item{evidence_sources}{the \code{evidence_sources} values
#'     \code{\link{get_network}()} filters on, as a list with
#'     \code{database} and \code{text_mined} elements, e.g.
#'     \code{\link{INDRA_DATABASE_SOURCES}}. STRING also has
#'     \code{predicted} (genomic context and co-expression channels).}
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
