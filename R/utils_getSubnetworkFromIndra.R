#' Validate input for MSstatsBioNet getSubnetworkFromIndra
#' @param input dataframe from MSstats groupComparison output
#' @param protein_level_data dataframe from MSstats dataProcess output
#' @param sources_filter sources filter
#' @param force_include_other character vector of identifiers to include in the
#' network.
#' @keywords internal
#' @noRd
.validateGetSubnetworkFromIndraInput <- function(input, protein_level_data, sources_filter, force_include_other) {
    if (!"HgncId" %in% colnames(input)) {
        stop("Invalid Input Error: Input must contain a column named 'HgncId'.")
    }
    num_proteins = length(unique(input$HgncId)) + 
        ifelse(!is.null(force_include_other), length(force_include_other), 0)
    if (num_proteins >= 400) {
        stop("Invalid Input Error: INDRA query must contain less than 400 proteins.  Consider lowering your p-value cutoff")
    }
    if (nrow(input) == 0) {
        stop("Invalid Input Error: Input must contain at least one protein after filtering.")
    }
    if (!is.null(protein_level_data)) {
        if(!all(c("Protein", "LogIntensities", "originalRUN") %in% colnames(protein_level_data))) {
            stop("protein_level_data must contain 'Protein', 'LogIntensities', and 'originalRUN' columns.")
        }
    }
    if (!is.null(sources_filter)) {
        if (!is.character(sources_filter)) {
            stop("sources_filter must be a character vector")
        }
    }
}

#' Call INDRA Cogex API and return response
#' @param hgncIds list of hgnc ids
#' @param force_include_other list of identifiers to include in the network
#' @return list of INDRA statements
#' @importFrom jsonlite toJSON
#' @importFrom httr POST add_headers content
#' @keywords internal
#' @noRd
.callIndraCogexApi <- function(hgncIds, force_include_other) {
    indraCogexUrl <-
        "https://discovery.indra.bio/api/indra_subnetwork_relations"

    hgncIds = unique(hgncIds)
    groundings <- lapply(hgncIds, function(x) list("HGNC", x))
    if (!is.null(force_include_other)) {
        groundings <- c(groundings, lapply(force_include_other, function(x) {
            parts <- unlist(strsplit(x, ":"))
            if (length(parts) != 2) {
                stop(paste0("Invalid identifier format: ", x, ". Expected format is 'namespace:identifier', e.g. 'HGNC:1234' or 'CHEBI:4911'."))
            }
            list(parts[1], parts[2])
        }))
    }
    groundings <- list(nodes = groundings)
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

#' @importFrom httr GET status_code content
#' @importFrom jsonlite fromJSON
.get_incorrect_curation_count <- function(stmt_hash, api_key) {
    stmt_hash_char <- as.character(stmt_hash)
    url <- paste0("https://db.indra.bio/curation/list/", stmt_hash_char, "?api_key=", api_key)

    tryCatch({
        response <- GET(url)
        if (status_code(response) == 200) {
            curations <- fromJSON(content(response, "text", encoding = "UTF-8"))
            if (length(curations) == 0) {
                return(0)
            }
            incorrect_curations <- curations[curations$tag != "correct", ]
            unique_incorrect <- length(unique(incorrect_curations$source_hash))
            
            return(unique_incorrect)
        } else {
            warning(paste("API request failed for hash", stmt_hash_char, 
                          "with status code", status_code(response)))
            return(0)
        }
    }, error = function(e) {
        warning(paste("Error processing hash", stmt_hash_char, ":", e$message))
        return(0)
    })
}

#' Call INDRA Cogex API and return response
#' @param res response from INDRA
#' @param interaction_types interaction types to filter by
#' @param evidence_count_cutoff number of evidence to filter on for each paper
#' @param sources_filter list of sources to filter by. Default is NULL, i.e. no filter
#' @param filter_by_curation logical, whether to filter out statements that
#' have been curated as incorrect in INDRA.  Default is FALSE.
#' @param api_key string of INDRA API key for accessing curated statements.
#' @return filtered list of INDRA statements
#' @importFrom jsonlite fromJSON
#' @keywords internal
#' @noRd
.filterIndraResponse <- function(res, interaction_types, evidence_count_cutoff, 
                                 sources_filter = NULL, filter_by_curation = FALSE, api_key = "") {
    if (!is.null(interaction_types)) {
        res = Filter(
            function(statement) statement$data$stmt_type %in% interaction_types, 
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
    if (filter_by_curation) {
        for (i in seq_along(res)) {
            stmt_json <- fromJSON(res[[i]]$data$stmt_json)
            stmt_hash <- stmt_json$matches_hash
            incorrect_count <- .get_incorrect_curation_count(stmt_hash, api_key)
            res[[i]]$data$evidence_count <- res[[i]]$data$evidence_count - incorrect_count
            # Todo: Also subtract source_counts accordingly if requested
            Sys.sleep(0.1)
        }
    }
    res = Filter(
        function(statement) statement$data$evidence_count >= evidence_count_cutoff, 
        res
    )
    return(res)
}

#' Filter groupComparison result input based on user-defined cutoffs
#' @param input groupComparison result
#' @param pvalueCutoff p-value cutoff
#' @param logfc_cutoff logFC cutoff
#' @param force_include_other list of identifiers to exempt from filtering
#' @return filtered groupComparison result
#' @keywords internal
#' @noRd
.filterGetSubnetworkFromIndraInput <- function(input, pvalueCutoff, logfc_cutoff, force_include_other) {
    input$Protein <- as.character(input$Protein)
    
    # Extract exempt proteins before any filtering
    exempt_proteins <- NULL
    if (!is.null(force_include_other)) {
        if (!is.character(force_include_other)) {
            stop("force_include_other must be a character vector")
        }
        if ("HgncId" %in% colnames(input) && any(grepl("^HGNC:", force_include_other))) {
            hgnc_ids_to_include <- gsub("^HGNC:", "", force_include_other[grepl("^HGNC:", force_include_other)])
            exempt_proteins <- input[input$HgncId %in% hgnc_ids_to_include, ]
        } else {
            exempt_proteins <- data.frame()
        }
    }
    
    # Apply standard filtering
    input <- input[!is.na(input$adj.pvalue),]
    if (!is.null(pvalueCutoff)) {
        input <- input[input$adj.pvalue < pvalueCutoff, ]
    }
    if (!is.null(logfc_cutoff)) {
        if (!is.numeric(logfc_cutoff) || length(logfc_cutoff) != 1 || logfc_cutoff < 0) {
            stop("logfc_cutoff must be a single positive numeric value")
        }
        input <- input[!is.na(input$log2FC) & abs(input$log2FC) > logfc_cutoff, ]
    }
    if ("issue" %in% colnames(input)) {
        input <- input[is.na(input$issue), ]
    }
    
    # Combine filtered data with exempt proteins and remove duplicates
    if (!is.null(exempt_proteins) && nrow(exempt_proteins) > 0) {
        combined_input <- rbind(exempt_proteins, input)
        # Remove duplicates based on Protein column, keeping first occurrence
        input <- combined_input[!duplicated(combined_input$Protein), ]
    }
    
    # Handle PTMs in Protein column
    input$Site = ifelse(grepl("_[A-Z][0-9]", input$Protein),
                        gsub("^_", "", 
                             gsub("^[^_]*_|_(?![A-Z][0-9])[^_]*", "", input$Protein, perl = TRUE)
                         ),
                        NA_character_
                )
    if ("GlobalProtein" %in% colnames(input)) {
        input$Protein = input$GlobalProtein
    }
    return(input)
}
#' Add additional metadata to an edge
#' @param edge object representation of an INDRA statement
#' @param input filtered groupComparison result
#' @return edge with additional metadata
#' @keywords internal
#' @noRd
.addAdditionalMetadataToIndraEdge <- function(edge, input) {
    edge$evidence_list <- paste(
        "https://db.indra.bio/statements/from_agents?subject=",
        edge$source_id, "@", edge$source_ns, "&object=",
        edge$target_id, "@", edge$target_ns, "&format=html",
        sep = ""
    )
    
    # Convert back to uniprot IDs
    matched_rows_source <- input[which(input$HgncId == edge$source_id), ]
    uniprot_ids_source <- unique(matched_rows_source$Protein)
    if (length(uniprot_ids_source) != 1) {
        edge$source_uniprot_id <- edge$source_name
    } else {
        edge$source_uniprot_id <- uniprot_ids_source
    }
    
    matched_rows_target <- input[which(input$HgncId == edge$target_id), ]
    uniprot_ids_target = unique(matched_rows_target$Protein)
    if (length(uniprot_ids_target) != 1) {
        edge$target_uniprot_id <- edge$target_name
    } else {
        edge$target_uniprot_id <- uniprot_ids_target
    }
    
    return(edge)
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
#' @param protein_level_data output of dataProcess
#' @importFrom r2r query keys
#' @importFrom jsonlite fromJSON
#' @return edge data.frame
#' @keywords internal
#' @noRd
.constructEdgesDataFrame <- function(res, input, protein_level_data) {
    res <- .collapseDuplicateEdgesIntoEdgeToMetadataMapping(res, input)
    edges <- data.frame(
        source = vapply(keys(res), function(x) {
            query(res, x)$source_uniprot_id
        }, ""),
        target = vapply(keys(res), function(x) {
            query(res, x)$target_uniprot_id
        }, ""),
        site = vapply(keys(res), function(x) {
            query(res, x)$site
        }, ""),
        interaction = vapply(keys(res), function(x) {
            query(res, x)$data$stmt_type
        }, ""),
        evidenceCount = vapply(keys(res), function(x) {
            query(res, x)$data$evidence_count
        }, 1),
        paperCount = vapply(keys(res), function(x) {
            query(res, x)$data$paper_count
        }, 1),
        evidenceLink = vapply(keys(res), function(x) {
            query(res, x)$evidence_list
        }, ""),
        sourceCounts = vapply(keys(res), function(x) {
            query(res, x)$data$source_counts
        }, ""),
        stmt_hash = vapply(keys(res), function(x) {
            stmt_json <- fromJSON(query(res, x)$data$stmt_json)
            stmt_json$matches_hash
        }, ""),
        stringsAsFactors = FALSE
    )
    # add correlation - maybe create a separate function
    if (!is.null(protein_level_data)) {
        protein_level_data <- protein_level_data[
            protein_level_data$Protein %in% edges$source | 
                protein_level_data$Protein %in% edges$target, ]
        correlations <- .getCorrelationMatrixFromProteinLevelData(protein_level_data)
        edges$correlation <- apply(edges, 1, function(edge) {
            if (edge["source"] %in% rownames(correlations) && edge["target"] %in% colnames(correlations)) {
                return(correlations[edge["source"], edge["target"]])
            } else {
                return(NA)
            }
        })
    }
    return(edges)
}

#' Construct nodes data.frame from groupComparison output
#' @param input filtered groupComparison result
#' @param edges edges data frame
#' @return nodes data.frame
#' @keywords internal
#' @noRd
.constructNodesDataFrame <- function(input, edges) {
    nodes = input[, c("Protein", "HgncName", "Site", "log2FC", "adj.pvalue")]
    colnames(nodes) = c("id", "hgncName", "Site", "logFC", "adj.pvalue")
    
    nodes = nodes[nodes$id %in% c(edges$source, edges$target), ]
    extra_force_include_other <- setdiff(unique(c(edges$source, edges$target)), nodes$id)
    if (length(extra_force_include_other) > 0) {
        extra_nodes <- data.frame(
            id = extra_force_include_other,
            hgncName = NA,
            Site = NA,
            logFC = 0,
            adj.pvalue = 1,
            stringsAsFactors = FALSE
        )
        nodes <- rbind(nodes, extra_nodes)
    }
    nodes$hgncName = ifelse(is.na(nodes$hgncName), nodes$id, nodes$hgncName)
    
    return(nodes)
}

#' Filter Edges Data Frame
#' @param edges response from INDRA
#' @param paper_count_cutoff cutoff for number of papers
#' @param correlation_cutoff if protein_level_abundance is not NULL, apply a 
#' cutoff for edges with correlation less than a specified cutoff.
#' @return filtered edges data frame
#' @keywords internal
#' @noRd
.filterEdgesDataFrame <- function(edges, 
                                  paper_count_cutoff,
                                  correlation_cutoff) {
    edges <- edges[which(edges$paperCount >= paper_count_cutoff), ]
    if ("correlation" %in% colnames(edges)) {
        edges <- edges[which(abs(edges$correlation) >= correlation_cutoff), ]
    }
    if (nrow(edges) == 0) {
        stop("No edges remain after applying filters. Consider relaxing filters")
    }
    return(edges)
}

#' Construct correlation matrix from MSstats
#' @param protein_level_data output of dataProcess
#' @importFrom tidyr pivot_wider
#' @importFrom stats cor
#' @return correlations matrix
#' @keywords internal
#' @noRd
.getCorrelationMatrixFromProteinLevelData <- function(protein_level_data) {
    Protein = LogIntensities = NULL
    wide_data <- pivot_wider(protein_level_data[,c("Protein", "LogIntensities", "originalRUN")], names_from = Protein, values_from = LogIntensities)
    wide_data <- wide_data[, -which(names(wide_data) == "originalRUN")]
    if (any(colSums(!is.na(wide_data)) == 0)) {
        warning("protein_level_data contains proteins with all missing values, unable to calculate correlations for those proteins.")
    }
    correlations <- cor(wide_data, use = "pairwise.complete.obs")
    return(correlations)
}
