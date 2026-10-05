# INDRA CoGEx HTTP calls and the parsers that turn their responses into
# the edge contract. Moved from utils_getSubnetworkFromIndra.R in Phase 2
# of the API refactor.

#' Build the list of (namespace, id) groundings for the INDRA Cogex query
#'
#' Splits each row's semicolon-joined \code{EntityNamespace} / \code{EntityId}
#' positionally, fans out each pair into its own grounding node, then appends
#' any \code{force_include_other} entries (parsed as \code{"namespace:id"}),
#' returning the unique set.
#' @param namespaces character vector aligned with \code{ids}
#' @param ids character vector aligned with \code{namespaces}
#' @param force_include_other optional character vector of
#'        \code{"namespace:id"} strings
#' @return list of two-element \code{list(namespace, id)} groundings
#' @keywords internal
#' @noRd
.buildCogexGroundings <- function(namespaces, ids, force_include_other = NULL) {
    if (length(namespaces) != length(ids)) {
        stop("EntityNamespace and EntityId must have the same length")
    }

    pairs <- list()
    for (i in seq_along(namespaces)) {
        ns_i <- strsplit(as.character(namespaces[i]), ";")[[1]]
        id_i <- strsplit(as.character(ids[i]),        ";")[[1]]
        if (length(ns_i) != length(id_i)) {
            stop("EntityNamespace and EntityId entries must be positionally aligned ",
                 "after splitting on ';' (mismatch at row ", i, ")")
        }
        for (k in seq_along(ns_i)) {
            pairs <- c(pairs, list(list(ns_i[k], id_i[k])))
        }
    }

    if (!is.null(force_include_other)) {
        for (x in force_include_other) {
            parts <- unlist(strsplit(x, ":"))
            if (length(parts) != 2) {
                stop(paste0("Invalid identifier format: ", x,
                            ". Expected format is 'namespace:identifier', e.g. 'HGNC:1234' or 'CHEBI:4911'."))
            }
            pairs <- c(pairs, list(list(parts[1], parts[2])))
        }
    }

    unique(pairs)
}

#' Call INDRA Cogex API and return response
#' @param namespaces character vector of entity namespaces (semicolon-joined
#'        for multi-grounded rows), aligned with \code{ids}
#' @param ids character vector of entity ids (semicolon-joined for
#'        multi-grounded rows), aligned with \code{namespaces}
#' @param force_include_other list of \code{"namespace:id"} identifiers to
#'        include in the network
#' @param cogex_url base URL of INDRA CoGEx
#' @return list of INDRA statements
#' @importFrom jsonlite toJSON
#' @importFrom httr POST add_headers content
#' @keywords internal
#' @noRd
.callIndraCogexApi <- function(namespaces, ids, force_include_other,
                               cogex_url = INDRA_API_URL) {
    indraCogexUrl <- file.path(cogex_url, "api/indra_subnetwork_relations")

    pairs <- .buildCogexGroundings(namespaces, ids, force_include_other)
    groundings <- list(nodes = pairs)
    groundings <- jsonlite::toJSON(groundings, auto_unbox = TRUE)

    res <- POST(
        indraCogexUrl,
        body = groundings,
        add_headers("Content-Type" = "application/json"),
        encode = "raw"
    )
    res <- content(res)
    return(res)
}

#' Call INDRA Cogex API and return response
#' @param res response from INDRA
#' @param statement_types interaction types to filter by
#' @param evidence_count_cutoff number of evidence to filter on for each paper
#' @param sources_filter list of sources to filter by. Default is NULL, i.e. no filter
#' @return filtered list of INDRA statements
#' @importFrom jsonlite fromJSON
#' @keywords internal
#' @noRd
.filterIndraResponse <- function(res, statement_types, evidence_count_cutoff, 
                                 sources_filter = NULL) {
    if (!is.null(statement_types)) {
        res = Filter(
            function(statement) statement$data$stmt_type %in% statement_types, 
            res)
    }
    if (!is.null(sources_filter)) {
        res = Filter(
            function(statement) {
                parsed <- tryCatch(fromJSON(statement$data$source_counts), error = function(e) NULL)
                if (is.null(parsed)) return(FALSE)
                return(any(names(parsed) %in% sources_filter))
            }, 
            res
        )
    }
    res = Filter(
        function(statement) statement$data$evidence_count >= evidence_count_cutoff, 
        res
    )
    return(res)
}

#' Row-level membership check for a (namespace, id) endpoint
#'
#' Splits each row's \code{EntityNamespace}/\code{EntityId} on \code{";"}
#' and returns \code{TRUE} for rows whose grounding list contains the
#' (\code{edge_ns}, \code{edge_id}) pair. Used to map an INDRA edge
#' endpoint back to the original \code{Protein} value.
#' @keywords internal
#' @noRd
.rowMatchesEndpoint <- function(input, edge_ns, edge_id) {
    row_ns <- strsplit(as.character(input$EntityNamespace), ";")
    row_id <- strsplit(as.character(input$EntityId),        ";")
    vapply(seq_along(row_ns), function(i) {
        rns <- row_ns[[i]]; rid <- row_id[[i]]
        if (length(rns) != length(rid) || length(rns) == 0) return(FALSE)
        any(rns == edge_ns & rid == edge_id)
    }, logical(1))
}

#' Add additional metadata to an edge
#' @param edge object representation of an INDRA statement
#' @param input filtered groupComparison result
#' @return edge with additional metadata
#' @keywords internal
#' @noRd
.addAdditionalMetadataToIndraEdge <- function(edge, input) {
    edge$evidence_url <- .indraStatementUrl(edge$statement_id)

    # Map the grounded INDRA endpoint back to the original Protein value.
    # Membership-test against each row's semicolon-split (namespace, id)
    # pairs, using INDRA's source_ns/target_ns for namespace-aware disambiguation.
    matched_rows_source <- input[.rowMatchesEndpoint(input, edge$source_ns, edge$source_id), ]
    node_ids_source <- unique(matched_rows_source$Protein)
    if (length(node_ids_source) != 1) {
        edge$source_node_id <- edge$source_name
    } else {
        edge$source_node_id <- node_ids_source
    }

    matched_rows_target <- input[.rowMatchesEndpoint(input, edge$target_ns, edge$target_id), ]
    node_ids_target <- unique(matched_rows_target$Protein)
    if (length(node_ids_target) != 1) {
        edge$target_node_id <- edge$target_name
    } else {
        edge$target_node_id <- node_ids_target
    }

    return(edge)
}

#' Build the INDRA DB page URL for one statement
#' @param statement_id INDRA statement hash, as character
#' @return character URL
#' @keywords internal
#' @noRd
.indraStatementUrl <- function(statement_id) {
    paste0("https://db.indra.bio/statements/from_hash/", statement_id,
           "?format=html")
}


#' Collapse duplicate INDRA statements into a mapping of edge to metadata
#' @param res INDRA response
#' @param input filtered groupComparison result
#' @importFrom jsonlite fromJSON
#' @importFrom r2r hashmap keys
#' @return processed edge to metadata mapping
#' @keywords internal
#' @noRd
.collapseDuplicateEdgesIntoEdgeToMetadataMapping <- function(res, input) {
    edgeToMetadataMapping <- hashmap()

    for (edge in res) {
        key <- paste(edge$source_id, edge$target_id, edge$data$stmt_type, sep = "_")
        json_object <- fromJSON(edge$data$stmt_json)
        # matches_hash is a JSON string, so it keeps full precision. The
        # numeric data$stmt_hash does not: hashes exceed 2^53.
        edge$statement_id <- as.character(json_object$matches_hash)
        if (!is.null(json_object$residue) && !is.null(json_object$position)) {
            edge$site = paste0(json_object$residue, json_object$position)
            key <- paste(key, edge$site, sep = "_")
        } else {
            edge$site = NA_character_
        }
        if (!key %in% keys(edgeToMetadataMapping) || 
            edge$data$evidence_count > edgeToMetadataMapping[[key]]$data$evidence_count) {
            edge <- .addAdditionalMetadataToIndraEdge(edge, input)
            edge$data$paper_count <- 1 # TODO: fix paper count
            edgeToMetadataMapping[[key]] <- edge
        }
    }

    return(edgeToMetadataMapping)
}

#' Construct edges data.frame from INDRA response
#' @param res INDRA response
#' @param input filtered groupComparison result
#' @importFrom r2r query keys
#' @importFrom jsonlite fromJSON
#' @return edge data.frame
#' @keywords internal
#' @noRd
.constructEdgesDataFrame <- function(res, input) {
    res <- .collapseDuplicateEdgesIntoEdgeToMetadataMapping(res, input)
    statements <- lapply(keys(res), function(x) query(res, x))
    interaction <- vapply(statements, function(x) x$data$stmt_type, "")
    edges <- data.frame(
        source = vapply(statements, function(x) x$source_node_id, ""),
        target = vapply(statements, function(x) x$target_node_id, ""),
        interaction = interaction,
        directed = !interaction %in% UNDIRECTED_STATEMENT_TYPES,
        site = vapply(statements, function(x) x$site, ""),
        confidence = vapply(statements, function(x) {
            if (is.null(x$data$belief)) NA_real_ else as.numeric(x$data$belief)
        }, 1),
        evidence_count = vapply(statements, function(x) {
            as.integer(x$data$evidence_count)
        }, 1L),
        evidence_url = vapply(statements, function(x) x$evidence_url, ""),
        statement_id = vapply(statements, function(x) x$statement_id, ""),
        backend_database = rep("INDRA", length(statements)),
        query_type = rep("subnetwork", length(statements)),
        evidence_sources = vapply(statements, function(x) {
            x$data$source_counts
        }, ""),
        paperCount = vapply(statements, function(x) x$data$paper_count, 1),
        stringsAsFactors = FALSE
    )
    return(edges)
}
