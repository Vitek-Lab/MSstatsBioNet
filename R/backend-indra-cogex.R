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

#' Add additional metadata to an edge
#' @param edge object representation of an INDRA statement
#' @param grounding_lookup the entity rows of each grounding, from
#' \code{.build_grounding_lookup()}
#' @return edge with additional metadata
#' @keywords internal
#' @noRd
.addAdditionalMetadataToIndraEdge <- function(edge, grounding_lookup) {
    edge$evidence_url <- .indraStatementUrl(edge$statement_id)
    # Map the statement's source and target back to the nodes of their
    # entity rows, matching namespace and identifier, so a multi-grounded row
    # is found by any of its groundings. A grounding can match several nodes.
    edge$source_node_ids <- .find_node_ids_for_grounding(
        grounding_lookup, edge$source_ns, edge$source_id, edge$source_name)
    edge$target_node_ids <- .find_node_ids_for_grounding(
        grounding_lookup, edge$target_ns, edge$target_id, edge$target_name)
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
#' @param grounding_lookup the entity rows of each grounding, from
#' \code{.build_grounding_lookup()}
#' @importFrom jsonlite fromJSON
#' @importFrom r2r hashmap keys
#' @return processed edge to metadata mapping
#' @keywords internal
#' @noRd
.collapseDuplicateEdgesIntoEdgeToMetadataMapping <- function(res, grounding_lookup) {
    edgeToMetadataMapping <- hashmap()

    for (edge in res) {
        # Identifiers are only unique within a namespace (HGNC:1234 is not
        # CHEBI:1234)
        key <- paste(edge$source_ns, edge$source_id, edge$target_ns,
                     edge$target_id, edge$data$stmt_type, sep = "_")
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
            edge <- .addAdditionalMetadataToIndraEdge(edge, grounding_lookup)
            edge$data$paper_count <- 1 # TODO: fix paper count
            edgeToMetadataMapping[[key]] <- edge
        }
    }

    return(edgeToMetadataMapping)
}

#' Construct edges data.frame from INDRA response
#' @param res INDRA response
#' @param grounding_lookup the entity rows of each grounding, from
#' \code{.build_grounding_lookup()}
#' @importFrom r2r query keys
#' @importFrom jsonlite fromJSON
#' @return edge data.frame
#' @keywords internal
#' @noRd
.constructEdgesDataFrame <- function(res, grounding_lookup) {
    res <- .collapseDuplicateEdgesIntoEdgeToMetadataMapping(res, grounding_lookup)
    statements <- lapply(keys(res), function(x) query(res, x))
    interaction <- vapply(statements, function(x) x$data$stmt_type, "")
    edges <- data.frame(
        source = character(length(statements)),
        target = character(length(statements)),
        interaction = interaction,
        directed = !interaction %in% UNDIRECTED_INTERACTION_TYPES,
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
    .fan_out_edges(edges, statements)
}

#' Give an edge to every pair of nodes its statement connects
#'
#' A statement whose source or target grounding matches several nodes
#' becomes one edge per (source node, target node) pair. The copies share
#' the statement's columns, including \code{statement_id} and
#' \code{evidence_count}. A statement from a grounding to itself (e.g. a
#' homodimer) gives each matching node a self-loop, and no edges between
#' those nodes, which the statement doesn't name.
#' @param edges edges data.frame, one row per statement, \code{source} and
#' \code{target} not yet filled in
#' @param statements the statements, with \code{source_node_ids} and
#' \code{target_node_ids} from \code{.addAdditionalMetadataToIndraEdge()}
#' @return edges data.frame, with the copies of a statement next to each
#' other
#' @keywords internal
#' @noRd
.fan_out_edges <- function(edges, statements) {
    node_pairs <- lapply(statements, function(statement) {
        source_ids <- statement$source_node_ids
        target_ids <- statement$target_node_ids
        same_grounding <- identical(statement$source_ns, statement$target_ns) &&
            identical(as.character(statement$source_id),
                      as.character(statement$target_id))
        if (same_grounding) {
            return(list(source = source_ids, target = source_ids))
        }
        list(source = rep(source_ids, each = length(target_ids)),
             target = rep(target_ids, times = length(source_ids)))
    })
    n_pairs <- vapply(node_pairs, function(x) length(x$source), integer(1))
    edges <- edges[rep(seq_len(nrow(edges)), n_pairs), , drop = FALSE]
    edges$source <- as.character(unlist(lapply(node_pairs, `[[`, "source")))
    edges$target <- as.character(unlist(lapply(node_pairs, `[[`, "target")))
    rownames(edges) <- NULL
    edges
}

# CoGEx ID-mapping and entity-property calls, used by convert_ids() and
# get_entity_properties(). Moved from utils_annotateProteinInfoFromIndra.R in
# Phase 3 of the API refactor.

#' Call API to get UniProt IDs from UniProt mnemonic IDs
#' @param uniprotMnemonicIds list of UniProt mnemonic ids
#' @param cogex_url base URL of INDRA CoGEx
#' @return list of UniProt IDs
#' @importFrom jsonlite toJSON
#' @importFrom httr POST add_headers content
#' @keywords internal
#' @noRd
.callGetUniprotIdsFromUniprotMnemonicIdsApi <- function(uniprotMnemonicIds, cogex_url = INDRA_API_URL) {

    if (!is.list(uniprotMnemonicIds)) {
        stop("Input must be a list.")
    }

    if (length(uniprotMnemonicIds) == 0) {
        stop("Input list must not be empty.")
    }

    tryCatch({
        # Attempt to convert all elements to character if not already character
        uniprotMnemonicIds <- lapply(uniprotMnemonicIds, function(x) {
            if (!is.character(x)) {
                as.character(x)
            } else {
                x
            }
        })

        # Check if conversion was successful
        if (any(!sapply(uniprotMnemonicIds, is.character))) {
            stop("All elements in the list must be character strings representing UniProt mnemonic IDs.")
        }
    }, error = function(e) {
        stop("An error occurred converting uniprot mnemonic IDs to character strings: ", e$message)
    })

    apiUrl <- file.path(cogex_url, "api/get_uniprot_ids_from_uniprot_mnemonic_ids")

    requestBody <- list(uniprot_mnemonic_ids = uniprotMnemonicIds)
    requestBody <- jsonlite::toJSON(requestBody, auto_unbox = TRUE)
    res <- tryCatch({
        response <- POST(
            apiUrl,
            body = requestBody,
            add_headers("Content-Type" = "application/json"),
            encode = "raw"
        )
        content(response)
    }, error = function(e) {
        message("Error in API call: ", e)
        NULL
    })
    return(res)
}

#' Call API to get HGNC IDs from UniProt IDs
#' @param uniprotIds list of UniProt IDs
#' @param cogex_url base URL of INDRA CoGEx
#' @return list of HGNC IDs
#' @importFrom jsonlite toJSON
#' @importFrom httr POST add_headers content
#' @keywords internal
#' @noRd
.callGetHgncIdsFromUniprotIdsApi <- function(uniprotIds, cogex_url = INDRA_API_URL) {

    if (!is.list(uniprotIds)) {
        stop("Input must be a list.")
    }

    if (any(!sapply(uniprotIds, is.character))) {
        stop("All elements in the list must be character strings representing UniProt IDs.")
    }

    if (length(uniprotIds) == 0) {
        stop("Input list must not be empty.")
    }

    apiUrl <- file.path(cogex_url, "api/get_hgnc_ids_from_uniprot_ids")

    requestBody <- list(uniprot_ids = uniprotIds)
    requestBody <- jsonlite::toJSON(requestBody, auto_unbox = TRUE)
    res <- tryCatch({
        response <- POST(
            apiUrl,
            body = requestBody,
            add_headers("Content-Type" = "application/json"),
            encode = "raw"
        )
        content(response)
    }, error = function(e) {
        message("Error in API call: ", e)
        NULL
    })
    return(res)
}

#' Call API to get HGNC names from HGNC IDs
#' @param hgncIds list of HGNC IDs
#' @param cogex_url base URL of INDRA CoGEx
#' @return list of HGNC names
#' @importFrom jsonlite toJSON
#' @importFrom httr POST add_headers content
#' @keywords internal
#' @noRd
.callGetHgncNamesFromHgncIdsApi <- function(hgncIds, cogex_url = INDRA_API_URL) {

    if (!is.list(hgncIds)) {
        stop("Input must be a list.")
    }

    if (any(!sapply(hgncIds, is.character))) {
        stop("All elements in the list must be character strings representing HGNC IDs.")
    }

    if (length(hgncIds) == 0) {
        stop("Input list must not be empty.")
    }

    apiUrl <- file.path(cogex_url, "api/get_hgnc_names_from_hgnc_ids")

    requestBody <- list(hgnc_ids = hgncIds)
    requestBody <- jsonlite::toJSON(requestBody, auto_unbox = TRUE)
    res <- tryCatch({
        response <- POST(
            apiUrl,
            body = requestBody,
            add_headers("Content-Type" = "application/json"),
            encode = "raw"
        )
        content(response)
    }, error = function(e) {
        message("Error in API call: ", e)
        NULL
    })
    return(res)
}

#' Call API to check if genes are kinases
#' @param genes list of gene names
#' @param cogex_url base URL of INDRA CoGEx
#' @return list indicating if genes are kinases
#' @importFrom jsonlite toJSON
#' @importFrom httr POST add_headers content
#' @keywords internal
#' @noRd
.callIsKinaseApi <- function(genes, cogex_url = INDRA_API_URL) {

    if (!is.list(genes)) {
        stop("Input must be a list.")
    }

    if (any(!sapply(genes, is.character))) {
        stop("All elements in the list must be character strings representing gene names.")
    }

    if (length(genes) == 0) {
        stop("Input list must not be empty.")
    }

    apiUrl <- file.path(cogex_url, "api/is_kinase")

    requestBody <- list(genes = genes)
    requestBody <- jsonlite::toJSON(requestBody, auto_unbox = TRUE)
    res <- tryCatch({
        response <- POST(
            apiUrl,
            body = requestBody,
            add_headers("Content-Type" = "application/json"),
            encode = "raw"
        )
        content(response)
    }, error = function(e) {
        message("Error in API call: ", e)
        NULL
    })
    return(res)
}

#' Call API to check if genes are phosphatases
#' @param genes list of gene names
#' @param cogex_url base URL of INDRA CoGEx
#' @return list indicating if genes are phosphatases
#' @importFrom jsonlite toJSON
#' @importFrom httr POST add_headers content
#' @keywords internal
#' @noRd
.callIsPhosphataseApi <- function(genes, cogex_url = INDRA_API_URL) {

    if (!is.list(genes)) {
        stop("Input must be a list.")
    }

    if (any(!sapply(genes, is.character))) {
        stop("All elements in the list must be character strings representing gene names.")
    }

    if (length(genes) == 0) {
        stop("Input list must not be empty.")
    }

    apiUrl <- file.path(cogex_url, "api/is_phosphatase")

    requestBody <- list(genes = genes)
    requestBody <- jsonlite::toJSON(requestBody, auto_unbox = TRUE)
    res <- tryCatch({
        response <- POST(
            apiUrl,
            body = requestBody,
            add_headers("Content-Type" = "application/json"),
            encode = "raw"
        )
        content(response)
    }, error = function(e) {
        message("Error in API call: ", e)
        NULL
    })
    return(res)
}

#' Call API to check if genes are transcription factors
#' @param genes list of gene names
#' @param cogex_url base URL of INDRA CoGEx
#' @return list indicating if genes are transcription factors
#' @importFrom jsonlite toJSON
#' @importFrom httr POST add_headers content
#' @keywords internal
#' @noRd
.callIsTranscriptionFactorApi <- function(genes, cogex_url = INDRA_API_URL) {

    if (!is.list(genes)) {
        stop("Input must be a list.")
    }

    if (any(!sapply(genes, is.character))) {
        stop("All elements in the list must be character strings representing gene names.")
    }

    if (length(genes) == 0) {
        stop("Input list must not be empty.")
    }

    apiUrl <- file.path(cogex_url, "api/is_transcription_factor")

    requestBody <- list(genes = genes)
    requestBody <- jsonlite::toJSON(requestBody, auto_unbox = TRUE)
    res <- tryCatch({
        response <- POST(
            apiUrl,
            body = requestBody,
            add_headers("Content-Type" = "application/json"),
            encode = "raw"
        )
        content(response)
    }, error = function(e) {
        message("Error in API call: ", e)
        NULL
    })
    return(res)
}


#' Query INDRA CoGEx for the evidence of statements
#'
#' Hashes that are missing from the response (unknown hashes, or a
#' batch whose request failed) are absent from the returned list.
#'
#' @param stmt_hashes Character vector of statement hash strings
#' @param cogex_url   base URL of INDRA CoGEx
#' @param batch_size  Number of hashes to request per API call
#' @param sleep       Seconds to pause between API calls
#' @return Named list: statement_id -> list of evidence objects. Empty list when
#'         no evidence could be retrieved.
#' @keywords internal
#' @noRd
#' @importFrom httr POST status_code content content_type_json
.query_indra_evidence <- function(stmt_hashes, cogex_url = INDRA_API_URL,
                                  batch_size = 100, sleep = 1) {
    url <- file.path(cogex_url, "api/get_evidences_for_stmt_hashes")

    stmt_hashes <- unique(as.character(stmt_hashes))
    if (length(stmt_hashes) == 0) return(list())

    batches   <- split(stmt_hashes, ceiling(seq_along(stmt_hashes) / batch_size))
    n_batches <- length(batches)
    results   <- list()

    cat(sprintf("Fetching evidence for %d hashes in %d batch(es) of up to %d...\n",
                length(stmt_hashes), n_batches, batch_size))

    for (i in seq_along(batches)) {
        batch <- batches[[i]]

        parsed <- tryCatch({
            response <- POST(
                url,
                body   = list(stmt_hashes = I(batch)),
                encode = "json",
                content_type_json()
            )

            if (status_code(response) != 200) {
                warning(sprintf("API returned status %d for stmt_hashes: %s",
                                status_code(response),
                                paste(batch, collapse = ", ")))
                NULL
            } else {
                content(response, as = "parsed")
            }
        }, error = function(e) {
            warning(sprintf("Error querying stmt_hashes %s: %s",
                            paste(batch, collapse = ", "), e$message))
            return(NULL)
        })

        cat(sprintf("Progress: %d/%d batches (%.1f%%)\n",
                    i, n_batches, (i / n_batches) * 100))

        if (!is.null(parsed)) {
            matched <- intersect(names(parsed), batch)
            results[matched] <- parsed[matched]
        }

        if (i < n_batches) Sys.sleep(sleep)
    }

    cat("Done fetching evidence!\n")
    results
}
