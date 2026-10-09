# STRING backend: protein-protein associations from the STRING database,
# through its REST API (https://string-db.org/help/api/). Phase 5d of the
# API refactor; the design is in TODO-MSBio-20261001_refactor_api.md, 4.8.

#' Address of STRING's current version, used to resolve the version
#' @keywords internal
#' @noRd
STRING_API_URL <- "https://string-db.org"

#' The namespace of STRING groundings, e.g. STRING:9606.ENSP00000269305
#' @keywords internal
#' @noRd
STRING_NAMESPACE <- "STRING"

#' STRING network types
#' @keywords internal
#' @noRd
STRING_NETWORK_TYPES <- c("physical", "functional", "regulatory")

#' Identifier systems the STRING backend converts, by entity type
#'
#' UniProt accessions only: get_string_ids is a text search, which maps a
#' mnemonic such as TP53_HUMAN to the wrong protein (HIPK4).
#' @keywords internal
#' @noRd
STRING_ID_CONVERSIONS <- list(
    protein  = "uniprot",
    ptm_site = "uniprot"
)

#' Pattern of a UniProt accession, with an optional isoform suffix
#'
#' From https://www.uniprot.org/help/accession_numbers.
#' @keywords internal
#' @noRd
UNIPROT_ACCESSION_PATTERN <- paste0(
    "^([OPQ][0-9][A-Z0-9]{3}[0-9]|",
    "[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})(-[0-9]+)?$")

#' Evidence channels of each STRING network type, and their response fields
#'
#' Names are the values of edges$evidence_sources and of the
#' evidence_sources argument of get_network().
#' @keywords internal
#' @noRd
STRING_CHANNELS <- list(
    physical = c(neighborhood = "nscore", fusion = "fscore",
                 cooccurence = "pscore", coexpression = "ascore",
                 experiments = "escore", database = "dscore",
                 textmining = "tscore"),
    regulatory = c(experiments = "experimental_score",
                   database = "database_score",
                   textmining = "textmining_score")
)
STRING_CHANNELS$functional <- STRING_CHANNELS$physical

#' STRING evidence channels, grouped as backend_capabilities() lists them
#'
#' "experiments" holds interactions imported from primary interaction
#' databases (BioGRID, IntAct, ...), and "database" from curated pathway
#' databases (KEGG, Reactome, ...), so both are database evidence. The
#' other channels are predictions from genomic context and co-expression.
#' @keywords internal
#' @noRd
STRING_CHANNEL_GROUPS <- list(
    database   = c("experiments", "database"),
    text_mined = "textmining",
    predicted  = c("neighborhood", "fusion", "cooccurence", "coexpression")
)

#' The interaction types of each STRING network type
#' @keywords internal
#' @noRd
STRING_INTERACTION_TYPES <- list(
    physical   = "Complex",
    functional = "Association",
    regulatory = c("Activation", "Inhibition", "Regulation")
)

#' Largest number of proteins one STRING query takes, by query type
#'
#' More than 10 proteins also need a species, which get_network() always
#' sends.
#' @keywords internal
#' @noRd
STRING_MAX_NODES <- c(subnetwork = 2000)

#' Create a STRING backend
#'
#' STRING is a database of protein-protein associations, scored by how
#' likely they are to be real. The backend queries STRING's REST API, so it
#' needs no other package. It grounds UniProt accessions to STRING proteins
#' and answers \code{\link{subnetwork_query}()}: the STRING edges among the
#' selected proteins. \code{backend_capabilities(string_backend())} lists
#' what it supports.
#'
#' STRING has three network types. \code{"physical"} (the default) has
#' edges between proteins that bind, as undirected \code{"Complex"} edges,
#' and is the closest to INDRA's mechanisms. \code{"functional"} adds
#' proteins that work together without binding, e.g. from co-expression
#' or text mining, as undirected \code{"Association"} edges, and is much
#' denser. \code{"regulatory"} has directed edges, as \code{"Activation"},
#' \code{"Inhibition"}, or \code{"Regulation"} (no sign).
#'
#' Each STRING version has its own address, and the backend sends every
#' request to the address of one version, so one backend never mixes
#' versions and \code{network$provenance} records the version. With
#' \code{version = NULL}, the current version is looked up on the first
#' request. STRING asks for no more than one request per second, so the
#' backend waits between requests.
#'
#' STRING has no evidence count: \code{evidence_count} is \code{NA}.
#' \code{evidence_sources} lists the evidence channels with a non-zero
#' score, e.g. \code{"database;experiments;textmining"}: \code{experiments}
#' (interactions imported from experimental databases such as BioGRID and
#' IntAct), \code{database} (curated pathway databases such as KEGG and
#' Reactome), \code{textmining} (co-mention in the literature), and, for
#' the physical and functional networks, the predictions
#' \code{neighborhood}, \code{fusion}, \code{cooccurence}, and
#' \code{coexpression}. A channel pools many databases, so STRING edges
#' can't be filtered by one database. \code{confidence} is STRING's
#' combined score of the channels. STRING has no PTM sites, evidence text, or
#' curations, so \code{\link{get_evidence}()} and
#' \code{\link{get_curations}()} give an error for it, and
#' \code{\link{filter_by_curation}()} keeps its edges as they are.
#'
#' STRING is licensed under the Creative Commons Attribution 4.0 license.
#' Cite STRING when using its networks; see
#' \url{https://string-db.org/cgi/about}.
#'
#' @param network_type \code{"physical"}, \code{"functional"}, or
#' \code{"regulatory"}
#' @param version STRING version, e.g. \code{"12.5"}. \code{NULL} uses the
#' current version.
#' @param caller_identity name sent to STRING with each request, so STRING
#' can see which tool its requests come from
#' @return a \code{StringBackend} object, to pass to
#' \code{\link{convert_ids}()} and \code{\link{get_network}()}
#' @seealso \code{\link{NetworkBackend-class}}, \code{\link{indra_backend}()}
#' @export
#' @examples
#' string <- string_backend()
#' backend_capabilities(string)$interaction_types
#' \donttest{
#' df <- data.frame(Protein = c("P04637", "Q00987", "P38936"),
#'                  adj.pvalue = c(0.01, 0.02, 0.03))
#' entities <- prepare_entities(df, entity_type = "protein",
#'                              id_type = "uniprot")
#' entities <- convert_ids(string, entities)
#' entities <- select_entities(entities, pvalue_cutoff = 0.05)
#' network <- get_network(string, entities, min_confidence = 0.7)
#' network$edges[, c("source", "target", "interaction", "confidence")]
#' }
string_backend <- function(network_type = c("physical", "functional",
                                            "regulatory"),
                           version = NULL,
                           caller_identity = "MSstatsBioNet") {
    network_type <- match.arg(network_type)
    new("StringBackend", network_type = network_type,
        version = if (is.null(version)) NA_character_ else
            as.character(version),
        caller_identity = caller_identity,
        host = new.env(parent = emptyenv()))
}

#' @rdname backend_capabilities
#' @export
setMethod("backend_capabilities", "StringBackend",
    function(backend) {
        channels <- names(STRING_CHANNELS[[backend@network_type]])
        list(query_types       = "subnetwork",
             id_conversions    = STRING_ID_CONVERSIONS,
             entity_properties = character(0),
             interaction_types =
                 STRING_INTERACTION_TYPES[[backend@network_type]],
             evidence_sources  = lapply(STRING_CHANNEL_GROUPS, intersect,
                                        channels),
             max_nodes         = STRING_MAX_NODES)
    })

# STRING identifier conversion. Sends the unique UniProt accessions of each
# organism to get_string_ids, without isoform suffixes, since STRING maps
# only canonical accessions. Rows of a protein group get the STRING protein
# of each member. Groundings of other backends (INDRA) are kept.
#' @rdname convert_ids
#' @export
setMethod("convert_ids", "StringBackend",
    function(backend, entities, ...) {
        .validate_entities(entities)
        .check_id_conversions(entities, backend, STRING_ID_CONVERSIONS)
        members <- .get_grounding_inputs(entities)
        .check_uniprot_accessions(unlist(members, use.names = FALSE))
        accessions <- lapply(members, function(row_members) {
            unique(sub("-[0-9]+$", "", row_members))
        })
        organisms <- .get_entity_organisms(entities)
        groundings <- .build_empty_groundings(nrow(entities))
        unmapped <- character(0)
        for (organism in unique(organisms)) {
            rows <- which(organisms == organism)
            queries <- unique(unlist(accessions[rows], use.names = FALSE))
            if (length(queries) == 0) {
                next
            }
            mapping <- .map_accessions_to_string(backend, queries, organism)
            unmapped <- c(unmapped, setdiff(queries, mapping$accession))
            for (i in rows) {
                matched <- mapping[mapping$accession %in% accessions[[i]], ,
                                   drop = FALSE]
                matched <- matched[order(match(matched$accession,
                                               accessions[[i]])), ,
                                   drop = FALSE]
                matched <- matched[!duplicated(matched$string_id), ,
                                   drop = FALSE]
                if (nrow(matched) == 0) {
                    next
                }
                groundings$namespace[i] <- .joinProteinGroup(
                    rep(STRING_NAMESPACE, nrow(matched)))
                groundings$entity_id[i] <- .joinProteinGroup(matched$string_id)
                groundings$entity_name[i] <- .joinProteinGroup(
                    matched$preferred_name)
            }
        }
        if (length(unmapped) > 0) {
            message(length(unmapped), " UniProt accession(s) not found in ",
                    "STRING: ", .list_values_for_message(unmapped), ".")
        }
        .replace_groundings(entities, seq_len(nrow(entities)), groundings,
                            .is_string_namespace)
    })

# STRING subnetwork query. Reads only the STRING groundings, sends the STRING
# proteins of the included_in_query rows to the network endpoint, which
# returns the edges among them, and normalizes the edges to the contract.
# Proteins from include_entities that match no row become latent nodes.
#' @rdname get_network
#' @export
setMethod("get_network", signature("StringBackend", "SubnetworkQuery"),
    function(backend, entities, query, interaction_types = NULL,
             min_evidence = 1, min_confidence = NULL,
             evidence_sources = NULL, include_entities = NULL, ...) {
        entities <- .keep_groundings(entities, .is_string_namespace)
        groundings <- .get_groundings_to_query(entities)
        .check_min_confidence(min_confidence)
        .validate_string_subnetwork_input(backend, groundings, min_evidence,
                                          evidence_sources, include_entities)
        organism <- .get_string_query_organism(entities)
        string_ids <- unique(c(
            groundings$entity_id,
            sub(paste0("^", STRING_NAMESPACE, ":"), "", include_entities)))
        message(.describe_subnetwork_question(
            entities, groundings, include_entities,
            backend_label = paste("STRING", backend@network_type)))
        host <- .resolve_string_host(backend)
        response <- .call_string_api(
            backend, "network",
            list(identifiers = paste(string_ids, collapse = "\r"),
                 species = organism,
                 network_type = backend@network_type,
                 required_score = .get_string_required_score(min_confidence)))
        retrieved_at <- .get_current_time()
        interactions <- .parse_string_network(response, backend@network_type)
        interactions <- .filter_string_interactions(
            interactions, interaction_types, evidence_sources,
            backend@network_type)
        grounding_lookup <- .build_grounding_lookup(entities)
        edges <- .build_string_edges(interactions, grounding_lookup,
                                     backend@network_type, host$address,
                                     organism)
        edges <- .filter_by_min_evidence(edges, min_evidence)
        edges <- .filter_by_min_confidence(edges, min_confidence)
        edges <- .filterEdgesDataFrame(edges)
        nodes <- .build_network_nodes(
            grounding_lookup, edges,
            .list_latent_string_nodes(interactions, grounding_lookup))
        nodes$entity_type[!nodes$measured &
                              nodes$namespace %in% STRING_NAMESPACE] <-
            "protein"
        provenance <- .build_provenance(
            backend_database = "STRING", query_type = "subnetwork",
            backend_url = host$address, organism = organism,
            parameters = list(interaction_types = interaction_types,
                              min_evidence = min_evidence,
                              min_confidence = min_confidence,
                              evidence_sources = evidence_sources,
                              include_entities = include_entities,
                              network_type = backend@network_type),
            backend_version = host$version,
            retrieved_at = retrieved_at)
        network <- list(nodes = nodes, edges = edges, provenance = provenance)
        validate_network(network)
        network
    })

#' Whether namespaces are STRING's
#' @param namespaces character vector
#' @return logical vector
#' @keywords internal
#' @noRd
.is_string_namespace <- function(namespaces) {
    namespaces %in% STRING_NAMESPACE
}

#' Stop unless every identifier is a UniProt accession
#'
#' Checked before any request: STRING's get_string_ids is a text search and
#' would map, e.g., a mnemonic in a column labelled "uniprot" to the wrong
#' protein.
#' @param identifiers character vector
#' @return \code{NULL}, invisibly; errors listing the others
#' @keywords internal
#' @noRd
.check_uniprot_accessions <- function(identifiers) {
    not_accessions <- unique(identifiers[!grepl(UNIPROT_ACCESSION_PATTERN,
                                                identifiers)])
    if (length(not_accessions) > 0) {
        stop("StringBackend takes UniProt accessions (e.g. \"P04637\"), ",
             "since STRING's ID search maps other text to the wrong ",
             "protein. Not accessions: ",
             .list_values_for_message(not_accessions), ".", call. = FALSE)
    }
    invisible(NULL)
}

#' Map UniProt accessions to STRING proteins
#'
#' STRING leaves out the identifiers it can't map, and also all but one of
#' several identifiers of the same protein, so the accessions missing from
#' the first response are sent once more on their own.
#' @param backend a \code{StringBackend}
#' @param accessions unique UniProt accessions, without isoform suffixes
#' @param organism NCBI taxon ID
#' @return data.frame with \code{accession}, \code{string_id},
#' \code{preferred_name}, for the accessions that map
#' @keywords internal
#' @noRd
.map_accessions_to_string <- function(backend, accessions, organism) {
    mapping <- .request_string_ids(backend, accessions, organism)
    missing <- setdiff(accessions, mapping$accession)
    if (length(missing) > 0 && length(missing) < length(accessions)) {
        mapping <- rbind(mapping,
                         .request_string_ids(backend, missing, organism))
    }
    mapping
}

#' Call get_string_ids for a set of accessions
#'
#' Each returned row is checked against the accession sent at its
#' \code{queryIndex} and against the organism, so a row STRING matched to
#' something else is dropped.
#' @inheritParams .map_accessions_to_string
#' @return data.frame like \code{.map_accessions_to_string()}'s
#' @keywords internal
#' @noRd
.request_string_ids <- function(backend, accessions, organism) {
    response <- .call_string_api(
        backend, "get_string_ids",
        list(identifiers = paste(accessions, collapse = "\r"),
             species = organism, echo_query = 1, limit = 1))
    mapping <- data.frame(accession = character(0), string_id = character(0),
                          preferred_name = character(0),
                          stringsAsFactors = FALSE)
    if (nrow(response) == 0) {
        return(mapping)
    }
    sent <- accessions[as.integer(response$queryIndex) + 1]
    matches_query <- !is.na(sent) & sent == response$queryItem &
        as.character(response$ncbiTaxonId) == organism
    rbind(mapping, data.frame(
        accession = sent[matches_query],
        string_id = as.character(response$stringId[matches_query]),
        preferred_name = as.character(response$preferredName[matches_query]),
        stringsAsFactors = FALSE))
}

#' Validate the input of the STRING subnetwork query
#' @param backend a \code{StringBackend}
#' @param groundings grounding table of the rows to query
#' @param min_evidence,evidence_sources,include_entities arguments of
#' \code{get_network()}
#' @keywords internal
#' @noRd
.validate_string_subnetwork_input <- function(backend, groundings,
                                              min_evidence, evidence_sources,
                                              include_entities) {
    if (nrow(groundings) == 0) {
        stop("No selected entity has a STRING grounding. Run ",
             "convert_ids(string_backend(), entities) first.", call. = FALSE)
    }
    n_proteins <- length(unique(groundings$entity_id)) +
        length(include_entities)
    if (n_proteins > STRING_MAX_NODES[["subnetwork"]]) {
        stop("A STRING query takes at most ",
             STRING_MAX_NODES[["subnetwork"]], " proteins; got ", n_proteins,
             ". Consider lowering your p-value cutoff.", call. = FALSE)
    }
    if (!is.numeric(min_evidence) || length(min_evidence) != 1 ||
        is.na(min_evidence)) {
        stop("min_evidence must be a single number.", call. = FALSE)
    }
    if (!is.null(evidence_sources)) {
        channels <- names(STRING_CHANNELS[[backend@network_type]])
        unknown <- setdiff(evidence_sources, channels)
        if (!is.character(evidence_sources) || length(unknown) > 0) {
            stop("evidence_sources for the STRING ", backend@network_type,
                 " network must be among: ", paste(channels, collapse = ", "),
                 ".", call. = FALSE)
        }
    }
    if (!is.null(include_entities)) {
        pattern <- paste0("^", STRING_NAMESPACE, ":.+")
        if (!is.character(include_entities) ||
            !all(grepl(pattern, include_entities))) {
            stop("include_entities for STRING must be STRING groundings, ",
                 "e.g. \"STRING:9606.ENSP00000269305\".", call. = FALSE)
        }
    }
    invisible(NULL)
}

#' The organism of a STRING query
#'
#' STRING takes one species per request.
#' @param entities entity table
#' @return NCBI taxon ID of the queried rows
#' @keywords internal
#' @noRd
.get_string_query_organism <- function(entities) {
    organisms <- unique(.get_entity_organisms(entities)[
        entities$included_in_query & !is.na(entities$entity_id)])
    if (length(organisms) != 1) {
        stop("A STRING query takes one organism; the selected entities ",
             "have: ", .list_values_for_message(organisms), ".",
             call. = FALSE)
    }
    organisms
}

#' The required_score of a STRING request
#'
#' STRING applies 0.4 when none is sent; \code{min_confidence = NULL} means
#' no cutoff, so 0 is sent.
#' @param min_confidence the argument of get_network()
#' @return integer in [0, 1000]
#' @keywords internal
#' @noRd
.get_string_required_score <- function(min_confidence) {
    if (is.null(min_confidence)) 0L else as.integer(round(1000 * min_confidence))
}

#' Parse a STRING network response
#'
#' @param response data.frame from the network endpoint
#' @param network_type the backend's network type
#' @return data.frame with one row per STRING interaction:
#' \code{source_id}, \code{target_id}, \code{source_name},
#' \code{target_name}, \code{interaction}, \code{interaction_raw},
#' \code{directed}, \code{sign}, \code{confidence}, \code{statement_id},
#' and one score column per evidence channel. Undirected interactions
#' listed in both orders are kept once.
#' @keywords internal
#' @noRd
.parse_string_network <- function(response, network_type) {
    channels <- STRING_CHANNELS[[network_type]]
    n <- nrow(response)
    if (network_type == "regulatory") {
        ids <- c("source_string_id", "target_string_id")
        names <- c("source_preferred_name", "target_preferred_name")
        score <- "combined_score"
    } else {
        ids <- c("stringId_A", "stringId_B")
        names <- c("preferredName_A", "preferredName_B")
        score <- "score"
    }
    if (n == 0) {
        response <- as.data.frame(
            stats::setNames(rep(list(character(0)), 4), c(ids, names)),
            stringsAsFactors = FALSE)
    }
    directed <- network_type == "regulatory"
    sign_raw <- if (directed) as.character(response$sign) else
        rep(NA_character_, n)
    interactions <- data.frame(
        source_id   = as.character(response[[ids[1]]]),
        target_id   = as.character(response[[ids[2]]]),
        source_name = as.character(response[[names[1]]]),
        target_name = as.character(response[[names[2]]]),
        interaction = .get_string_interaction(network_type, sign_raw),
        interaction_raw = if (directed) {
            paste0("regulatory:", sign_raw, recycle0 = TRUE)
        } else {
            rep(network_type, n)
        },
        directed    = rep(directed, n),
        sign        = ifelse(sign_raw %in% "pos", 1L,
                             ifelse(sign_raw %in% "neg", -1L, NA_integer_)),
        confidence  = as.numeric(response[[score]]),
        stringsAsFactors = FALSE
    )
    for (channel in names(channels)) {
        values <- response[[channels[[channel]]]]
        interactions[[channel]] <- if (is.null(values)) rep(0, n) else
            as.numeric(values)
    }
    pair <- if (directed) {
        paste(interactions$source_id, interactions$target_id, sep = "|")
    } else {
        paste(pmin(interactions$source_id, interactions$target_id),
              pmax(interactions$source_id, interactions$target_id), sep = "|")
    }
    interactions$statement_id <- paste0("string:", network_type, ":", pair,
                                        recycle0 = TRUE)
    interactions[!duplicated(interactions$statement_id), , drop = FALSE]
}

#' The interaction type of STRING edges
#' @param network_type the backend's network type
#' @param sign_raw STRING's \code{sign} of each edge (regulatory only)
#' @return character vector
#' @keywords internal
#' @noRd
.get_string_interaction <- function(network_type, sign_raw) {
    switch(network_type,
        physical   = rep("Complex", length(sign_raw)),
        functional = rep("Association", length(sign_raw)),
        regulatory = ifelse(sign_raw %in% "pos", "Activation",
                            ifelse(sign_raw %in% "neg", "Inhibition",
                                   "Regulation")))
}

#' Filter STRING interactions by type and evidence channel
#' @param interactions from \code{.parse_string_network()}
#' @param interaction_types,evidence_sources arguments of get_network()
#' @param network_type the backend's network type
#' @return \code{interactions}, filtered. With \code{evidence_sources}, an
#' interaction is kept when one of the named channels has a non-zero score.
#' @keywords internal
#' @noRd
.filter_string_interactions <- function(interactions, interaction_types,
                                        evidence_sources, network_type) {
    if (!is.null(interaction_types)) {
        interactions <- interactions[
            interactions$interaction %in% interaction_types, , drop = FALSE]
    }
    if (!is.null(evidence_sources)) {
        scores <- as.matrix(interactions[, evidence_sources, drop = FALSE])
        interactions <- interactions[rowSums(scores > 0) > 0, , drop = FALSE]
    }
    interactions
}

#' Build the edges of a STRING network
#'
#' Each interaction becomes one edge per (source node, target node) pair,
#' as for INDRA, so rows that share a STRING protein (two isoforms of one
#' protein) each get its edges.
#' @param interactions from \code{.parse_string_network()}
#' @param grounding_lookup from \code{.build_grounding_lookup()}
#' @param network_type the backend's network type
#' @param address the STRING version's address, for \code{evidence_url}
#' @param organism NCBI taxon ID, for \code{evidence_url}
#' @return edges data.frame that meets the contract
#' @keywords internal
#' @noRd
.build_string_edges <- function(interactions, grounding_lookup,
                                network_type, address, organism) {
    n <- nrow(interactions)
    channels <- names(STRING_CHANNELS[[network_type]])
    evidence_sources <- vapply(seq_len(n), function(i) {
        scores <- unlist(interactions[i, channels, drop = TRUE])
        .join_evidence_sources(names(scores)[scores > 0])
    }, character(1))
    edges <- data.frame(
        source           = character(n),
        target           = character(n),
        interaction      = interactions$interaction,
        directed         = interactions$directed,
        site             = rep(NA_character_, n),
        confidence       = interactions$confidence,
        evidence_count   = rep(NA_integer_, n),
        evidence_url     = paste0(address, "/cgi/network?identifiers=",
                                  interactions$source_id, "%0d",
                                  interactions$target_id, "&species=",
                                  organism, recycle0 = TRUE),
        statement_id     = interactions$statement_id,
        backend_database = rep("STRING", n),
        query_type       = rep("subnetwork", n),
        evidence_sources = evidence_sources,
        interaction_raw  = interactions$interaction_raw,
        sign             = interactions$sign,
        stringsAsFactors = FALSE
    )
    statements <- lapply(seq_len(n), function(i) {
        list(source_ns = STRING_NAMESPACE,
             source_id = interactions$source_id[i],
             target_ns = STRING_NAMESPACE,
             target_id = interactions$target_id[i],
             source_node_ids = .find_node_ids_for_grounding(
                 grounding_lookup, STRING_NAMESPACE,
                 interactions$source_id[i], interactions$source_name[i]),
             target_node_ids = .find_node_ids_for_grounding(
                 grounding_lookup, STRING_NAMESPACE,
                 interactions$target_id[i], interactions$target_name[i]))
    })
    .add_edge_per_node_pair(edges, statements)
}

#' List the STRING proteins of interactions that match no entity row
#' @param interactions from \code{.parse_string_network()}
#' @param grounding_lookup from \code{.build_grounding_lookup()}
#' @return data.frame like \code{.list_latent_backend_nodes()}'s, with
#' STRING's preferred name as \code{id}
#' @keywords internal
#' @noRd
.list_latent_string_nodes <- function(interactions, grounding_lookup) {
    backend_nodes <- data.frame(
        namespace   = rep(STRING_NAMESPACE, 2 * nrow(interactions)),
        entity_id   = c(interactions$source_id, interactions$target_id),
        entity_name = c(interactions$source_name, interactions$target_name),
        stringsAsFactors = FALSE)
    matches_no_row <- vapply(backend_nodes$entity_id, function(entity_id) {
        length(.find_entity_rows_for_grounding(
            grounding_lookup, STRING_NAMESPACE, entity_id)) == 0
    }, logical(1), USE.NAMES = FALSE)
    backend_nodes <- backend_nodes[matches_no_row, , drop = FALSE]
    backend_nodes$id <- backend_nodes$entity_name
    backend_nodes[!duplicated(backend_nodes$id), , drop = FALSE]
}

#' Resolve the version and address a STRING backend sends requests to
#'
#' Resolved once per backend object and kept in its \code{host}
#' environment: a given version has a fixed address
#' (\code{https://version-12-5.string-db.org}), and \code{version = NULL}
#' asks STRING's version endpoint for the current one.
#' @param backend a \code{StringBackend}
#' @return list with \code{version} and \code{address}
#' @keywords internal
#' @noRd
.resolve_string_host <- function(backend) {
    host <- backend@host
    if (is.null(host$address)) {
        if (is.na(backend@version)) {
            response <- .send_string_request(
                paste0(STRING_API_URL, "/api/json/version"),
                list(caller_identity = backend@caller_identity))
            host$version <- as.character(response$string_version[1])
            host$address <- as.character(response$stable_address[1])
        } else {
            host$version <- backend@version
            host$address <- paste0("https://version-",
                                   gsub(".", "-", backend@version,
                                        fixed = TRUE),
                                   ".string-db.org")
        }
    }
    list(version = host$version, address = host$address)
}

#' Call an endpoint of the backend's STRING version
#' @param backend a \code{StringBackend}
#' @param endpoint e.g. \code{"network"}
#' @param parameters named list of request parameters
#' @return data.frame of the JSON response, with no rows when it is empty
#' @keywords internal
#' @noRd
.call_string_api <- function(backend, endpoint, parameters) {
    address <- .resolve_string_host(backend)$address
    .send_string_request(
        paste0(address, "/api/json/", endpoint),
        c(parameters, caller_identity = backend@caller_identity))
}

#' When the last STRING request was sent, for the one-per-second limit
#' @keywords internal
#' @noRd
.string_request_times <- new.env(parent = emptyenv())

#' Send a request to STRING, at most one per second
#'
#' STRING asks for no more than one request per second and no parallel
#' requests.
#' @param url endpoint URL
#' @param body named list of form parameters
#' @return data.frame of the JSON response, with no rows when it is empty
#' @keywords internal
#' @noRd
.send_string_request <- function(url, body) {
    last_request <- .string_request_times$last
    if (!is.null(last_request)) {
        elapsed <- as.numeric(difftime(.get_string_clock(), last_request,
                                       units = "secs"))
        if (elapsed < 1) {
            .wait_for_string(1 - elapsed)
        }
    }
    on.exit(.string_request_times$last <- .get_string_clock())
    response <- .post_string_request(url, body)
    if (length(response) == 0) {
        return(data.frame())
    }
    as.data.frame(response, stringsAsFactors = FALSE)
}

#' The clock for the STRING request limit; its own function for tests
#' @keywords internal
#' @noRd
.get_string_clock <- function() {
    Sys.time()
}

#' Wait before the next STRING request; its own function for tests
#' @param seconds how long
#' @keywords internal
#' @noRd
.wait_for_string <- function(seconds) {
    Sys.sleep(seconds)
}

#' POST a form to STRING and parse the JSON response
#' @param url endpoint URL
#' @param body named list of form parameters
#' @return the parsed JSON, a data.frame or an empty list
#' @importFrom httr POST content http_error status_code
#' @keywords internal
#' @noRd
.post_string_request <- function(url, body) {
    response <- httr::POST(url, body = body, encode = "form")
    if (httr::http_error(response)) {
        stop("STRING request to ", url, " failed with HTTP status ",
             httr::status_code(response), ".", call. = FALSE)
    }
    jsonlite::fromJSON(httr::content(response, as = "text",
                                     encoding = "UTF-8"))
}
