#' Base URL of INDRA CoGEx
#' @keywords internal
#' @noRd
INDRA_API_URL <- "https://discovery.indra.bio"

#' Base URL of Gilda, INDRA's grounding service
#' @keywords internal
#' @noRd
GILDA_API_URL <- "https://grounding.indra.bio"

#' Base URL of the INDRA database, which holds the curations
#' @keywords internal
#' @noRd
INDRA_DB_URL <- "https://db.indra.bio"

#' Create an INDRA backend
#'
#' INDRA is a knowledge graph of mechanisms (activations, phosphorylations,
#' complexes, ...) assembled from the literature and curated databases. The
#' backend queries INDRA CoGEx for networks and grounds gene symbols and
#' chemical names with Gilda, INDRA's grounding service. Curations, which
#' mark evidence as correct or incorrect, come from the INDRA database.
#' \code{backend_capabilities(indra_backend())} lists what it supports.
#'
#' This function includes third-party software components that are licensed
#' under the BSD 2-Clause License. Include the third-party licensing
#' agreements if redistributing this package or results based on it. See
#' the LICENSE file for details.
#'
#' @param cogex_url base URL of INDRA CoGEx
#' @param grounding_url base URL of Gilda
#' @param curation_url base URL of the INDRA database, which holds the
#' curations
#' @return an \code{IndraBackend} object, to pass to
#' \code{\link{convert_ids}()}, \code{\link{get_entity_properties}()},
#' \code{\link{get_network}()}, \code{\link{get_evidence}()}, and
#' \code{\link{get_curations}()}
#' @seealso \code{\link{NetworkBackend-class}}
#' @importFrom methods new
#' @export
#' @examples
#' indra <- indra_backend()
#' backend_capabilities(indra)$query_types
indra_backend <- function(cogex_url = INDRA_API_URL,
                          grounding_url = GILDA_API_URL,
                          curation_url = INDRA_DB_URL) {
    new("IndraBackend", cogex_url = cogex_url, grounding_url = grounding_url,
        curation_url = curation_url)
}

# INDRA subnetwork query. Sends the groundings of the included_in_query rows
# to CoGEx indra_subnetwork_relations, filters the statements, and
# normalizes them to the edge contract. The source and target of each
# statement are matched back to the rows of entities by grounding, so nodes
# carry their statistics. Backend nodes that match no row (from
# include_entities) become latent nodes.
#' @rdname get_network
#' @export
setMethod("get_network", signature("IndraBackend", "SubnetworkQuery"),
    function(backend, entities, query, interaction_types = NULL,
             min_evidence = 1, min_confidence = NULL,
             evidence_sources = NULL, include_entities = NULL, ...) {
        groundings <- .get_groundings_to_query(entities)
        .validateIndraSubnetworkInput(groundings, evidence_sources,
                                      include_entities)
        .check_min_confidence(min_confidence)
        message(.describe_subnetwork_question(entities, groundings,
                                              include_entities))
        statements <- .callIndraCogexApi(groundings$namespace,
                                         groundings$entity_id,
                                         include_entities, backend@cogex_url)
        statements <- .filterIndraResponse(statements, interaction_types,
                                           min_evidence, evidence_sources)
        grounding_lookup <- .build_grounding_lookup(entities)
        edges <- .constructEdgesDataFrame(statements, grounding_lookup)
        edges <- .filter_by_min_confidence(edges, min_confidence)
        edges <- .filterEdgesDataFrame(edges)
        nodes <- .build_network_nodes(
            grounding_lookup, edges,
            .list_backend_nodes(statements, grounding_lookup))
        network <- list(nodes = nodes, edges = edges)
        validate_network(network)
        network
    })

#' Largest number of groundings one INDRA query takes, by query type
#'
#' CoGEx \code{indra_subnetwork_relations} takes fewer than 400, so at most
#' 399.
#' @keywords internal
#' @noRd
INDRA_MAX_NODES <- c(subnetwork = 399)

#' @rdname backend_capabilities
#' @export
setMethod("backend_capabilities", "IndraBackend",
    function(backend) {
        list(query_types       = "subnetwork",
             id_conversions    = INDRA_ID_CONVERSIONS,
             entity_properties = names(INDRA_ENTITY_PROPERTIES),
             interaction_types = INTERACTION_TYPES,
             max_nodes         = INDRA_MAX_NODES)
    })

# INDRA evidence. Sends the unique statement hashes to CoGEx
# get_evidences_for_stmt_hashes and copies each evidence sentence onto every
# edge row with that hash. Network Search edges carry the same hashes, so
# this covers every INDRA query type.
#' @rdname get_evidence
#' @export
setMethod("get_evidence", "IndraBackend",
    function(backend, edges, ...) {
        .check_evidence_edge_columns(edges)
        statement_ids <- unique(edges$statement_id)
        cat(sprintf("Processing %d unique statement hashes...\n",
                    length(statement_ids)))
        evidence_by_statement <- .query_indra_evidence(statement_ids,
                                                       backend@cogex_url)
        .build_indra_evidence_table(edges, evidence_by_statement)
    })

# INDRA curations. Asks the INDRA database for the curations of each unique
# statement hash, one request per hash, and counts the evidence curated as
# anything other than correct.
#' @rdname get_curations
#' @export
setMethod("get_curations", "IndraBackend",
    function(backend, edges, ...) {
        if (!"statement_id" %in% names(edges)) {
            stop("Missing required columns: statement_id")
        }
        statement_ids <- unique(as.character(edges$statement_id))
        incorrect_counts <- integer(length(statement_ids))
        for (i in seq_along(statement_ids)) {
            if (i > 1) Sys.sleep(0.1)
            incorrect_counts[i] <- as.integer(.get_incorrect_curation_count(
                statement_ids[i], backend@curation_url))
        }
        data.frame(statement_id = statement_ids,
                   incorrect_count = incorrect_counts,
                   stringsAsFactors = FALSE)
    })

#' Edge columns that get_evidence() copies to its output
#' @keywords internal
#' @noRd
EVIDENCE_EDGE_COLUMNS <- c("source", "target", "interaction", "site",
                           "evidence_url", "statement_id")

#' Check that edges have the columns get_evidence() needs
#' @param edges edges data.frame
#' @return \code{NULL}, invisibly; errors naming the missing columns
#' @keywords internal
#' @noRd
.check_evidence_edge_columns <- function(edges) {
    missing_cols <- setdiff(EVIDENCE_EDGE_COLUMNS, names(edges))
    if (length(missing_cols) > 0) {
        stop(sprintf("Missing required columns: %s",
                     paste(missing_cols, collapse = ", ")))
    }
    invisible(NULL)
}

#' Build the evidence table from the CoGEx evidence of each statement
#' @param edges edges data.frame with \code{EVIDENCE_EDGE_COLUMNS}
#' @param evidence_by_statement named list from \code{.query_indra_evidence()}
#' @return data.frame with \code{EVIDENCE_EDGE_COLUMNS}, \code{text}, and
#'   \code{pmid}; one row per (edge, evidence with text) pair
#' @keywords internal
#' @noRd
.build_indra_evidence_table <- function(edges, evidence_by_statement) {
    results_list <- list()
    result_count <- 0

    for (statement_id in unique(edges$statement_id)) {
        evidence_list <- evidence_by_statement[[as.character(statement_id)]]
        if (is.null(evidence_list) || length(evidence_list) == 0) next

        matching_indices <- which(edges$statement_id == statement_id)

        for (evidence in evidence_list) {
            if (!is.null(evidence[["text"]]) && nchar(evidence[["text"]]) > 0) {
                for (idx in matching_indices) {
                    result_count <- result_count + 1
                    results_list[[result_count]] <- data.frame(
                        source       = edges$source[idx],
                        target       = edges$target[idx],
                        interaction  = edges$interaction[idx],
                        site         = edges$site[idx],
                        evidence_url = edges$evidence_url[idx],
                        statement_id = edges$statement_id[idx],
                        text         = evidence[["text"]],
                        pmid         = if (is.null(evidence[["pmid"]])) "" else evidence[["pmid"]],
                        stringsAsFactors = FALSE
                    )
                }
            }
        }
    }

    if (result_count == 0) {
        warning("No evidence text found for any statement hash")
        return(.build_empty_evidence_table())
    }

    results_df <- do.call(rbind, results_list)
    cat(sprintf("\nComplete! Found %d evidence text entries.\n", nrow(results_df)))
    results_df
}

#' Evidence table with no rows
#' @return data.frame with the get_evidence() columns and no rows
#' @keywords internal
#' @noRd
.build_empty_evidence_table <- function() {
    data.frame(
        source = character(), target = character(), interaction = character(),
        site = character(), evidence_url = character(), statement_id = character(),
        text = character(), pmid = character(), stringsAsFactors = FALSE
    )
}

#' The question a subnetwork query asks, with the entity counts
#'
#' Printed by \code{get_network()}, e.g. "INDRA subnetwork: how are 42
#' selected proteins connected to each other, with no other nodes added?".
#' @param entities entity table
#' @param groundings groundings of the rows to query, from
#' \code{.get_groundings_to_query()}
#' @param include_entities groundings added to the query
#' @return a single string
#' @keywords internal
#' @noRd
.describe_subnetwork_question <- function(entities, groundings,
                                          include_entities) {
    queried_types <- entities$entity_type[entities$id %in% groundings$id]
    selected <- .describe_entity_count(queried_types, "selected")
    if (length(include_entities) > 0) {
        selected <- paste0(selected, " and ", length(include_entities),
                           " added ",
                           if (length(include_entities) == 1) "entity"
                           else "entities")
    }
    paste0("INDRA subnetwork: how are ", selected, " connected to each ",
           "other, with no other nodes added?")
}

#' Names of entity types in messages, singular and plural
#' @keywords internal
#' @noRd
ENTITY_TYPE_NAMES <- list(
    protein    = c("protein", "proteins"),
    gene       = c("gene", "genes"),
    transcript = c("transcript", "transcripts"),
    ptm_site   = c("PTM site", "PTM sites"),
    metabolite = c("metabolite", "metabolites"),
    lipid      = c("lipid", "lipids"),
    drug       = c("drug", "drugs"),
    complex    = c("complex", "complexes"),
    family     = c("family", "families"),
    other      = c("entity", "entities")
)

#' Count entities by type in words, e.g. "42 selected proteins"
#' @param entity_types entity type of each entity
#' @param adjective word before the type, e.g. "selected"
#' @return a single string. Several types are counted as "entities".
#' @keywords internal
#' @noRd
.describe_entity_count <- function(entity_types, adjective) {
    types <- unique(entity_types)
    names <- if (length(types) == 1) ENTITY_TYPE_NAMES[[types]] else
        ENTITY_TYPE_NAMES$other
    noun <- if (length(entity_types) == 1) names[1] else names[2]
    paste(format(length(entity_types), big.mark = ","), adjective, noun)
}

#' Groundings of the entities to query
#'
#' Rows with \code{included_in_query} but no grounding can't be sent to a
#' backend, so they are left out with a message.
#' @param entities entity table
#' @return the \code{build_grounding_table()} of the \code{included_in_query}
#' rows
#' @keywords internal
#' @noRd
.get_groundings_to_query <- function(entities) {
    .validate_entities(entities)
    queried_rows <- entities[entities$included_in_query, , drop = FALSE]
    ungrounded <- is.na(queried_rows$namespace) |
        is.na(queried_rows$entity_id)
    if (any(ungrounded)) {
        message("Dropping ", sum(ungrounded),
                " row(s) with no entity grounding (NA entity_id).")
    }
    build_grounding_table(queried_rows[!ungrounded, , drop = FALSE])
}

#' Validate the input of the INDRA subnetwork query
#' @param groundings grounding table of the rows to query, from
#' \code{.get_groundings_to_query()}
#' @param evidence_sources evidence sources filter
#' @param include_entities character vector of \code{"namespace:identifier"}
#' groundings to add to the query
#' @keywords internal
#' @noRd
.validateIndraSubnetworkInput <- function(groundings, evidence_sources,
                                          include_entities) {
    unique_groundings <- unique(paste(groundings$namespace,
                                      groundings$entity_id, sep = ":"))
    num_proteins <- length(unique_groundings) + length(include_entities)
    if (num_proteins > INDRA_MAX_NODES[["subnetwork"]]) {
        stop("Invalid Input Error: INDRA query must contain less than 400 proteins.  Consider lowering your p-value cutoff")
    }
    if (nrow(groundings) == 0) {
        stop("Invalid Input Error: Input must contain at least one protein after filtering.")
    }
    if (!is.null(evidence_sources)) {
        if (!is.character(evidence_sources)) {
            stop("evidence_sources must be a character vector")
        }
    }
}

#' Entity types of the namespaces INDRA returns
#'
#' Used for latent nodes, which have no entity row to take a type from.
#' Namespaces not listed are \code{"other"}.
#' @keywords internal
#' @noRd
INDRA_NAMESPACE_ENTITY_TYPES <- c(
    HGNC      = "protein",
    UP        = "protein",
    UPPRO     = "protein",
    FPLX      = "family",
    CHEBI     = "metabolite",
    HMDB      = "metabolite",
    PUBCHEM   = "metabolite",
    LIPIDMAPS = "lipid",
    CHEMBL    = "drug",
    DRUGBANK  = "drug"
)

#' The entity type of a namespace
#' @param namespaces character vector of namespaces
#' @return character vector of entity types, \code{"other"} for unknown or
#' \code{NA} namespaces
#' @keywords internal
#' @noRd
.get_namespace_entity_types <- function(namespaces) {
    types <- unname(INDRA_NAMESPACE_ENTITY_TYPES[as.character(namespaces)])
    types[is.na(types)] <- "other"
    types
}

#' The node ID of each entity row
#'
#' A PTM site is drawn on its parent protein's node, so \code{ptm_site} rows
#' get their \code{parent_id}. Other rows keep their \code{id}.
#' @param entities entity table
#' @return character vector, one per row
#' @keywords internal
#' @noRd
.get_node_ids <- function(entities) {
    node_ids <- entities$id
    if ("parent_id" %in% colnames(entities)) {
        use_parent <- entities$entity_type == "ptm_site" &
            !is.na(entities$parent_id)
        node_ids[use_parent] <- entities$parent_id[use_parent]
    }
    node_ids
}

#' Look up the entity rows of each grounding
#'
#' Built once per query, so that matching each statement's source and target
#' to their entity rows doesn't split the \code{";"}-joined groundings again.
#' @param entities entity table
#' @return list with \code{entities}, \code{node_ids} (one per entity row,
#' from \code{.get_node_ids()}), and \code{rows_by_grounding}: for each
#' \code{"namespace:identifier"} grounding, the numbers of the entity rows
#' that have it
#' @keywords internal
#' @noRd
.build_grounding_lookup <- function(entities) {
    grounding_table <- build_grounding_table(entities)
    grounding_keys <- paste(grounding_table$namespace,
                            grounding_table$entity_id, sep = ":")
    list(entities = entities,
         node_ids = .get_node_ids(entities),
         rows_by_grounding = split(match(grounding_table$id, entities$id),
                                   grounding_keys))
}

#' Find the entity rows that have a grounding
#'
#' Prefers rows in the query. Rows outside it are matched only when no
#' queried row has the grounding, so that an unselected row with the same
#' grounding (an isoform, say) doesn't make the match ambiguous.
#' @param grounding_lookup from \code{.build_grounding_lookup()}
#' @param namespace,entity_id the grounding
#' @return integer vector of entity row numbers, empty when none match
#' @keywords internal
#' @noRd
.find_entity_rows_for_grounding <- function(grounding_lookup, namespace,
                                            entity_id) {
    rows <- grounding_lookup$rows_by_grounding[[
        paste(namespace, entity_id, sep = ":")]]
    if (is.null(rows)) {
        return(integer(0))
    }
    queried_rows <- rows[grounding_lookup$entities$included_in_query[rows]]
    if (length(queried_rows) > 0) queried_rows else rows
}

#' Find the node ID for a grounding returned by a backend
#' @param grounding_lookup from \code{.build_grounding_lookup()}
#' @param namespace,entity_id the grounding
#' @param backend_name the backend's name for the entity, e.g. INDRA's
#' gene symbol
#' @return the node ID of the matching entity rows, or \code{backend_name}
#' when none match or they belong to several nodes
#' @keywords internal
#' @noRd
.find_node_id_for_grounding <- function(grounding_lookup, namespace,
                                        entity_id, backend_name) {
    matching_rows <- .find_entity_rows_for_grounding(grounding_lookup,
                                                     namespace, entity_id)
    node_ids <- unique(grounding_lookup$node_ids[matching_rows])
    if (length(node_ids) == 1) node_ids else backend_name
}

#' List the backend nodes of INDRA statements
#' @param statements filtered INDRA response
#' @param grounding_lookup from \code{.build_grounding_lookup()}
#' @return data.frame with one row per distinct node \code{id}: the
#' sources and targets of the statements, with INDRA's \code{namespace},
#' \code{entity_id}, and \code{entity_name}
#' @keywords internal
#' @noRd
.list_backend_nodes <- function(statements, grounding_lookup) {
    backend_nodes <- lapply(statements, function(statement) {
        data.frame(
            namespace   = c(statement$source_ns, statement$target_ns),
            entity_id   = as.character(c(statement$source_id,
                                         statement$target_id)),
            entity_name = c(statement$source_name, statement$target_name),
            stringsAsFactors = FALSE)
    })
    backend_nodes <- do.call(rbind, c(list(.build_empty_groundings(0)),
                                      backend_nodes))
    backend_nodes$id <- vapply(seq_len(nrow(backend_nodes)), function(i) {
        .find_node_id_for_grounding(grounding_lookup,
                                    backend_nodes$namespace[i],
                                    backend_nodes$entity_id[i],
                                    backend_nodes$entity_name[i])
    }, character(1))
    backend_nodes[!duplicated(backend_nodes$id), , drop = FALSE]
}

#' Column order of the nodes that get_network() returns
#' @keywords internal
#' @noRd
NODE_COLUMN_ORDER <- c("id", "entity_type", "entity_name", "namespace",
                       "entity_id", "measured", "included_in_query",
                       "node_role", "site", "has_measured_sites", "logFC",
                       "adj.pvalue")

#' Build the nodes of a subnetwork from the entity table
#'
#' Nodes that an edge connects get one row per entity row: the rows in the
#' query, or, for a node that is there only through
#' \code{include_entities}, all of its rows. These are \code{measured}, with
#' their statistics. The remaining sources and targets are built by
#' \code{.build_latent_and_ambiguous_nodes()}.
#' @param grounding_lookup from \code{.build_grounding_lookup()}
#' @param edges edges data.frame
#' @param backend_nodes from \code{.list_backend_nodes()}
#' @return nodes data.frame with the columns in \code{NODE_COLUMN_ORDER}
#' @keywords internal
#' @noRd
.build_network_nodes <- function(grounding_lookup, edges, backend_nodes) {
    entities <- grounding_lookup$entities
    node_ids <- grounding_lookup$node_ids
    node_ids_in_edges <- unique(c(edges$source, edges$target))
    node_ids_in_query <- unique(node_ids[entities$included_in_query])
    keep_row <- node_ids %in% node_ids_in_edges &
        (entities$included_in_query | !node_ids %in% node_ids_in_query)
    entity_rows <- entities[keep_row, , drop = FALSE]
    user_added <- .get_column_or_default(entity_rows, "user_added",
                                         FALSE) %in% TRUE
    nodes_in_input <- data.frame(
        id                = node_ids[keep_row],
        entity_type       = entity_rows$entity_type,
        entity_name       = entity_rows$entity_name,
        namespace         = entity_rows$namespace,
        entity_id         = entity_rows$entity_id,
        measured          = rep(TRUE, nrow(entity_rows)),
        included_in_query = rep(TRUE, nrow(entity_rows)),
        node_role         = ifelse(entity_rows$included_in_query & !user_added,
                                   "passed_cutoffs", "user_added"),
        site              = .get_column_or_default(entity_rows, "site",
                                                   NA_character_),
        logFC             = .get_column_or_default(entity_rows, "logFC",
                                                   NA_real_),
        adj.pvalue        = .get_column_or_default(entity_rows, "adj.pvalue",
                                                   NA_real_),
        stringsAsFactors  = FALSE
    )
    unmatched <- backend_nodes$id %in% setdiff(node_ids_in_edges, node_ids)
    unmatched_backend_nodes <- backend_nodes[unmatched, , drop = FALSE]
    nodes <- rbind(nodes_in_input,
                   .build_latent_and_ambiguous_nodes(grounding_lookup,
                                                     unmatched_backend_nodes))
    nodes$has_measured_sites <- nodes$id %in%
        node_ids[entities$entity_type == "ptm_site"]
    nodes$entity_name <- ifelse(is.na(nodes$entity_name), nodes$id,
                                nodes$entity_name)
    nodes <- nodes[, NODE_COLUMN_ORDER]
    rownames(nodes) <- NULL
    nodes
}

#' Build nodes for backend nodes that are not the node of one entity
#'
#' A backend node that matches no entity row is latent: \code{measured =
#' FALSE}, \code{NA} statistics, and an \code{entity_type} from its
#' namespace. In a subnetwork query it can only come from
#' \code{include_entities}, so it is \code{"user_added"}. A backend node
#' that matches the rows of several nodes is ambiguous: it is in the input, but
#' it isn't known which row's statistics apply, so they are \code{NA}, and a
#' message names it.
#' @param grounding_lookup from \code{.build_grounding_lookup()}
#' @param backend_nodes rows of \code{.list_backend_nodes()}
#' @return nodes data.frame without \code{has_measured_sites}
#' @keywords internal
#' @noRd
.build_latent_and_ambiguous_nodes <- function(grounding_lookup, backend_nodes) {
    entities <- grounding_lookup$entities
    matching_rows <- lapply(seq_len(nrow(backend_nodes)), function(i) {
        .find_entity_rows_for_grounding(grounding_lookup,
                                        backend_nodes$namespace[i],
                                        backend_nodes$entity_id[i])
    })
    matches_several_nodes <- lengths(matching_rows) > 0
    if (any(matches_several_nodes)) {
        ambiguous_ids <- backend_nodes$id[matches_several_nodes]
        message(length(ambiguous_ids), " node(s) from the backend ",
                "match entities of several nodes and are shown under the ",
                "backend's name, without statistics: ",
                .list_values_for_message(ambiguous_ids), ".")
    }
    entity_types <- .get_namespace_entity_types(backend_nodes$namespace)
    node_roles <- rep("user_added", nrow(backend_nodes))
    user_added <- .get_column_or_default(entities, "user_added",
                                         FALSE) %in% TRUE
    for (i in which(matches_several_nodes)) {
        rows <- matching_rows[[i]]
        matching_types <- unique(entities$entity_type[rows])
        if (length(matching_types) == 1) {
            entity_types[i] <- matching_types
        }
        if (any(entities$included_in_query[rows] & !user_added[rows])) {
            node_roles[i] <- "passed_cutoffs"
        }
    }
    n_backend_nodes <- nrow(backend_nodes)
    data.frame(
        id                = backend_nodes$id,
        entity_type       = entity_types,
        entity_name       = backend_nodes$entity_name,
        namespace         = backend_nodes$namespace,
        entity_id         = backend_nodes$entity_id,
        measured          = matches_several_nodes,
        included_in_query = rep(TRUE, n_backend_nodes),
        node_role         = node_roles,
        site              = rep(NA_character_, n_backend_nodes),
        logFC             = rep(NA_real_, n_backend_nodes),
        adj.pvalue        = rep(NA_real_, n_backend_nodes),
        stringsAsFactors  = FALSE
    )
}

#' A column of a table, or a default value when the table has no such column
#' @keywords internal
#' @noRd
.get_column_or_default <- function(table, column, default) {
    if (column %in% colnames(table)) {
        table[[column]]
    } else {
        rep(default, nrow(table))
    }
}

#' Split a protein group into its member identifiers
#'
#' A \code{Protein} value may name a protein group -- several identifiers
#' for the same quantified analyte joined by \code{";"}. Splits on
#' \code{";"}, trims surrounding whitespace and drops empty members, so a
#' plain single identifier comes back as a length-one vector.
#'
#' @param x A length-one character value, possibly \code{NA}.
#' @return A character vector of member identifiers, empty when the input
#'         holds none.
#' @keywords internal
#' @noRd
.splitProteinGroup <- function(x) {
        if (length(x) == 0 || is.na(x)) {
                return(character(0))
        }
        members <- trimws(unlist(strsplit(as.character(x), ";", fixed = TRUE),
                                 use.names = FALSE))
        return(members[nzchar(members)])
}

#' Join protein group members back into a single value
#'
#' @param members A character vector of member identifiers.
#' @return The members joined by \code{";"}, or \code{NA} when empty.
#' @keywords internal
#' @noRd
.joinProteinGroup <- function(members) {
        if (length(members) == 0) {
                return(NA_character_)
        }
        return(paste(members, collapse = ";"))
}

#' Identifier systems the INDRA backend converts, by entity type
#'
#' Proteins (and the parent proteins of PTM sites) are grounded to HGNC
#' through CoGEx or Gilda. Chemical names are grounded through Gilda to
#' whichever namespace it returns (CHEBI, PUBCHEM, CHEMBL, ...).
#' @keywords internal
#' @noRd
INDRA_ID_CONVERSIONS <- list(
    protein    = c("uniprot", "uniprot_mnemonic", "hgnc_symbol"),
    ptm_site   = c("uniprot", "uniprot_mnemonic", "hgnc_symbol"),
    metabolite = "chemical_name",
    lipid      = "chemical_name",
    drug       = "chemical_name"
)

#' Entity properties of the INDRA backend
#'
#' \code{api} names the CoGEx call, which takes a list of HGNC gene symbols.
#' A PTM site gets the properties of its parent protein.
#' @keywords internal
#' @noRd
INDRA_ENTITY_PROPERTIES <- list(
    is_transcription_factor = list(api = ".callIsTranscriptionFactorApi",
                                   entity_types = c("protein", "ptm_site")),
    is_kinase               = list(api = ".callIsKinaseApi",
                                   entity_types = c("protein", "ptm_site")),
    is_phosphatase          = list(api = ".callIsPhosphataseApi",
                                   entity_types = c("protein", "ptm_site"))
)

# INDRA identifier conversion. Groups rows by id_type and makes one batch of
# calls per group: uniprot through CoGEx's UniProt-to-HGNC mapping,
# uniprot_mnemonic through CoGEx's mnemonic-to-UniProt mapping first,
# hgnc_symbol through Gilda restricted to HGNC and the row's organism, and
# chemical_name through Gilda with no namespace restriction. Each
# ";"-separated member of an identifier (a protein group) is grounded on its
# own, and the groundings are pooled onto the row.
#' @rdname convert_ids
#' @export
setMethod("convert_ids", "IndraBackend",
    function(backend, entities, ...) {
        .validate_entities(entities)
        .check_indra_id_conversions(entities)
        members <- .get_grounding_inputs(entities)
        for (id_type in unique(entities$id_type)) {
            rows <- which(entities$id_type == id_type)
            groundings <- switch(id_type,
                uniprot = .ground_uniprot_with_cogex(
                    members[rows], backend@cogex_url),
                uniprot_mnemonic = .ground_uniprot_with_cogex(
                    .map_uniprot_mnemonics(members[rows], backend@cogex_url),
                    backend@cogex_url),
                hgnc_symbol = .ground_text_with_gilda(
                    members[rows], backend@grounding_url, keep_only = "HGNC",
                    organisms = as.list(unique(.get_entity_organisms(entities)[rows]))),
                chemical_name = .ground_text_with_gilda(
                    members[rows], backend@grounding_url))
            entities$namespace[rows] <- groundings$namespace
            entities$entity_id[rows] <- groundings$entity_id
            entities$entity_name[rows] <- groundings$entity_name
        }
        entities
    })

# INDRA entity properties. Supports the properties in
# INDRA_ENTITY_PROPERTIES. They are looked up by gene symbol, so only rows
# with a single HGNC grounding get values; rows with several groundings (a
# protein group, or an ambiguous name) are NA.
#' @rdname get_entity_properties
#' @export
setMethod("get_entity_properties", "IndraBackend",
    function(backend, entities, properties = NULL, ...) {
        .validate_entities(entities)
        properties <- .resolve_entity_properties(
            backend, properties, names(INDRA_ENTITY_PROPERTIES))
        # A ";"-joined namespace (several groundings) never equals "HGNC"
        single_hgnc <- entities$namespace %in% "HGNC" &
            !is.na(entities$entity_name)
        for (property in properties) {
            spec <- INDRA_ENTITY_PROPERTIES[[property]]
            queried <- single_hgnc & entities$entity_type %in% spec$entity_types
            values <- rep(NA, nrow(entities))
            genes <- unique(entities$entity_name[queried])
            if (length(genes) > 0) {
                # Looked up by name at call time, so tests can mock the call
                call_api <- get(spec$api, mode = "function")
                response <- call_api(as.list(genes), backend@cogex_url)
                for (gene in names(response)) {
                    if (!is.null(response[[gene]])) {
                        values[queried & entities$entity_name == gene] <-
                            response[[gene]]
                    }
                }
            }
            entities[[property]] <- values
        }
        entities
    })

#' Stop when the INDRA backend can't convert some rows
#'
#' Checks each row's (\code{entity_type}, \code{id_type}) pair against
#' \code{INDRA_ID_CONVERSIONS}, and that proteins are human: CoGEx maps
#' UniProt IDs to HGNC, and Gilda is restricted to HGNC, which covers human
#' genes only.
#' @param entities entity table
#' @keywords internal
#' @noRd
.check_indra_id_conversions <- function(entities) {
    pairs <- unique(entities[, c("entity_type", "id_type")])
    supported <- vapply(seq_len(nrow(pairs)), function(i) {
        pairs$id_type[i] %in% INDRA_ID_CONVERSIONS[[pairs$entity_type[i]]]
    }, logical(1))
    if (any(!supported)) {
        unsupported <- paste(pairs$entity_type[!supported],
                             pairs$id_type[!supported], sep = " / ")
        allowed <- unlist(lapply(names(INDRA_ID_CONVERSIONS), function(type) {
            paste(type, INDRA_ID_CONVERSIONS[[type]], sep = " / ")
        }))
        stop("IndraBackend can't convert entity_type / id_type: ",
             .list_values_for_message(unsupported), ". Supported: ",
             paste(allowed, collapse = ", "), ".", call. = FALSE)
    }
    is_protein <- entities$entity_type %in% c("protein", "ptm_site")
    organisms <- unique(.get_entity_organisms(entities)[is_protein])
    non_human <- setdiff(organisms, "9606")
    if (length(non_human) > 0) {
        stop("IndraBackend grounds proteins to HGNC, which covers human ",
             "(organism \"9606\") only. Got organism: ",
             .list_values_for_message(non_human), ".", call. = FALSE)
    }
}

#' The organism of each entity row
#' @param entities entity table
#' @return character vector, \code{"9606"} where the table has no
#' \code{organism} column
#' @keywords internal
#' @noRd
.get_entity_organisms <- function(entities) {
    if ("organism" %in% colnames(entities)) {
        as.character(entities$organism)
    } else {
        rep("9606", nrow(entities))
    }
}

#' The identifiers to ground for each entity row
#' @param entities entity table
#' @return list with one character vector per row: the members of
#' \code{id}, or of \code{parent_id} for \code{ptm_site} rows
#' @keywords internal
#' @noRd
.get_grounding_inputs <- function(entities) {
    inputs <- entities$id
    if ("parent_id" %in% colnames(entities)) {
        use_parent <- entities$entity_type == "ptm_site" &
            !is.na(entities$parent_id)
        inputs[use_parent] <- entities$parent_id[use_parent]
    }
    lapply(inputs, .splitProteinGroup)
}

#' Map UniProt mnemonics to UniProt IDs through CoGEx
#' @param members list of character vectors of mnemonics, one per row
#' @param cogex_url base URL of INDRA CoGEx
#' @return list of character vectors of UniProt IDs, one per row. Members
#' that don't map are dropped.
#' @keywords internal
#' @noRd
.map_uniprot_mnemonics <- function(members, cogex_url) {
    mnemonics <- unique(unlist(members, use.names = FALSE))
    if (length(mnemonics) == 0) {
        return(members)
    }
    mapping <- .callGetUniprotIdsFromUniprotMnemonicIdsApi(as.list(mnemonics),
                                                           cogex_url)
    lapply(members, function(row_members) {
        unique(as.character(unlist(mapping[row_members], use.names = FALSE)))
    })
}

#' Ground UniProt IDs to HGNC through CoGEx
#' @param members list of character vectors of UniProt IDs, one per row
#' @param cogex_url base URL of INDRA CoGEx
#' @return data.frame with one row per element of \code{members}:
#' \code{namespace}, \code{entity_id}, \code{entity_name}, \code{";"}-joined
#' and \code{NA} for rows that don't ground. \code{entity_name} has
#' \code{"NA"} where one of several name lookups failed, and is \code{NA}
#' when all of them did.
#' @keywords internal
#' @noRd
.ground_uniprot_with_cogex <- function(members, cogex_url) {
    groundings <- .build_empty_groundings(length(members))
    uniprot_ids <- unique(unlist(members, use.names = FALSE))
    if (length(uniprot_ids) == 0) {
        return(groundings)
    }
    hgnc_mapping <- .callGetHgncIdsFromUniprotIdsApi(as.list(uniprot_ids),
                                                     cogex_url)
    hgnc_ids <- unique(as.character(unlist(hgnc_mapping, use.names = FALSE)))
    name_mapping <- list()
    if (length(hgnc_ids) > 0) {
        name_response <- .callGetHgncNamesFromHgncIdsApi(as.list(hgnc_ids),
                                                         cogex_url)
        if (!is.null(name_response)) {
            name_mapping <- name_response
        }
    }
    for (i in seq_along(members)) {
        entity_ids <- unique(as.character(
            unlist(hgnc_mapping[members[[i]]], use.names = FALSE)))
        if (length(entity_ids) == 0) {
            next
        }
        entity_names <- vapply(name_mapping[entity_ids], function(name) {
            if (is.null(name)) NA_character_ else as.character(name)[1]
        }, character(1), USE.NAMES = FALSE)
        groundings$namespace[i] <-
            .joinProteinGroup(rep("HGNC", length(entity_ids)))
        groundings$entity_id[i] <- .joinProteinGroup(entity_ids)
        if (!all(is.na(entity_names))) {
            groundings$entity_name[i] <- .joinProteinGroup(entity_names)
        }
    }
    groundings
}

#' Ground names through Gilda
#' @param members list of character vectors of names, one per row
#' @param grounding_url base URL of Gilda
#' @param keep_only namespace to keep, or \code{NULL} to keep all
#' @param organisms list of NCBI taxon IDs to restrict Gilda to, or
#' \code{NULL}
#' @return data.frame like \code{.ground_uniprot_with_cogex()}'s. Every
#' candidate Gilda returns is kept, in its ranking order, without duplicate
#' (namespace, identifier) pairs.
#' @keywords internal
#' @noRd
.ground_text_with_gilda <- function(members, grounding_url, keep_only = NULL,
                                    organisms = NULL) {
    groundings <- .build_empty_groundings(length(members))
    texts <- unique(unlist(members, use.names = FALSE))
    if (length(texts) == 0) {
        return(groundings)
    }
    grounding_map <- .callGroundEntitiesFromGildaApi(
        as.list(texts), keep_only = keep_only, organisms = organisms,
        grounding_url = grounding_url)
    if (is.null(grounding_map)) {
        return(groundings)
    }
    for (i in seq_along(members)) {
        namespaces <- character(0)
        entity_ids <- character(0)
        entity_names <- character(0)
        for (candidates in grounding_map[members[[i]]]) {
            if (is.null(candidates)) {
                next
            }
            stopifnot(length(candidates$ns) == length(candidates$id),
                      length(candidates$ns) == length(candidates$name))
            namespaces <- c(namespaces, as.character(candidates$ns))
            entity_ids <- c(entity_ids, as.character(candidates$id))
            entity_names <- c(entity_names, as.character(candidates$name))
        }
        keep <- !duplicated(paste(namespaces, entity_ids, sep = ":"))
        if (!any(keep)) {
            next
        }
        groundings$namespace[i] <- .joinProteinGroup(namespaces[keep])
        groundings$entity_id[i] <- .joinProteinGroup(entity_ids[keep])
        groundings$entity_name[i] <- .joinProteinGroup(entity_names[keep])
    }
    groundings
}

#' Ungrounded grounding columns for n rows
#' @keywords internal
#' @noRd
.build_empty_groundings <- function(n) {
    data.frame(namespace = rep(NA_character_, n),
               entity_id = rep(NA_character_, n),
               entity_name = rep(NA_character_, n),
               stringsAsFactors = FALSE)
}

#' Check the properties asked of get_entity_properties()
#' @param backend the backend, for the message
#' @param properties requested properties, or \code{NULL} for all
#' @param supported the properties the backend supports
#' @return the properties to add
#' @keywords internal
#' @noRd
.resolve_entity_properties <- function(backend, properties, supported) {
    if (is.null(properties)) {
        return(supported)
    }
    if (!is.character(properties) || anyNA(properties)) {
        stop("properties must be a character vector.", call. = FALSE)
    }
    unknown <- setdiff(properties, supported)
    if (length(unknown) > 0) {
        stop(class(backend), " does not support these entity properties: ",
             .list_values_for_message(unknown), ". Supported: ",
             paste(supported, collapse = ", "), ".", call. = FALSE)
    }
    unique(properties)
}
