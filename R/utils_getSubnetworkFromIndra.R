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
    num_proteins = nrow(input) + ifelse(!is.null(force_include_other), length(force_include_other), 0)
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

#' Call INDRA Cogex API and return response
#' @param res response from INDRA
#' @param interaction_types interaction types to filter by
#' @param evidence_count_cutoff number of evidence to filter on for each paper
#' @param sources_filter list of sources to filter by. Default is NULL, i.e. no filter
#' @return filtered list of INDRA statements
#' @importFrom jsonlite fromJSON
#' @keywords internal
#' @noRd
.filterIndraResponse <- function(res, interaction_types, evidence_count_cutoff, sources_filter = NULL) {
    if (!is.null(interaction_types)) {
        res = Filter(
            function(statement) statement$data$stmt_type %in% interaction_types, 
            res)
    }
    res = Filter(
        function(statement) statement$data$evidence_count >= evidence_count_cutoff, 
        res
    )
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
    return(res)
}

#' Filter groupComparison result input based on user-defined cutoffs
#' @param input groupComparison result
#' @param pvalueCutoff p-value cutoff
#' @param logfc_cutoff logFC cutoff
#' @param force_include_proteins list of proteins to exempt from filtering
#' @return filtered groupComparison result
#' @keywords internal
#' @noRd
.filterGetSubnetworkFromIndraInput <- function(input, pvalueCutoff, logfc_cutoff, force_include_proteins) {
    # Extract exempt proteins before any filtering
    exempt_proteins <- NULL
    if (!is.null(force_include_proteins)) {
        if (!is.character(force_include_proteins)) {
            stop("force_include_proteins must be a character vector")
        }
        missing_prots <- setdiff(force_include_proteins, input$Protein)
        if (length(missing_prots) > 0) {
            warning("force_include_proteins not found: ", paste(missing_prots, collapse = ", "))
        }
        exempt_proteins <- input[input$Protein %in% force_include_proteins,]
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
    
    input$Protein <- as.character(input$Protein)
    # TODO: Make PTM processing more robust
    if ("GlobalProtein" %in% colnames(input)) {
        input$Site <- as.character(input$Protein)
        input$Protein <- as.character(input$GlobalProtein)
        # Summarize rows by protein: mean adj.pvalue and log2FC, keep other columns
        if (any(duplicated(input$Protein))) {
            # Get unique proteins
            unique_proteins <- unique(input$Protein)
            
            # Initialize result data frame
            result_list <- list()
            
            for (i in seq_along(unique_proteins)) {
                protein <- unique_proteins[i]
                protein_rows <- input[input$Protein == protein, , drop = FALSE]
                
                if (nrow(protein_rows) == 1) {
                    # Single row, keep as is
                    result_list[[i]] <- protein_rows
                } else {
                    # Multiple rows, summarize
                    summarized_row <- protein_rows[1, , drop = FALSE]  # Start with first row
                    
                    # Calculate mean for adj.pvalue if it exists
                    if ("adj.pvalue" %in% names(protein_rows) && is.numeric(protein_rows$adj.pvalue)) {
                        summarized_row$adj.pvalue <- mean(protein_rows$adj.pvalue, na.rm = TRUE)
                    }
                    
                    # Calculate mean for log2FC if it exists
                    if ("log2FC" %in% names(protein_rows) && is.numeric(protein_rows$log2FC)) {
                        summarized_row$log2FC <- mean(protein_rows$log2FC, na.rm = TRUE)
                    }
                    
                    result_list[[i]] <- summarized_row
                }
            }
            
            # Combine all results
            input <- do.call(rbind, result_list)
            rownames(input) <- NULL
        }
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
    if (nrow(matched_rows_source) != 1) {
        edge$source_uniprot_id <- edge$source_name
    } else {
        edge$source_uniprot_id <- matched_rows_source$Protein
    }
    
    matched_rows_target <- input[which(input$HgncId == edge$target_id), ]
    if (nrow(matched_rows_target) != 1) {
        edge$target_uniprot_id <- edge$target_name
    } else {
        edge$target_uniprot_id <- matched_rows_target$Protein
    }
    
    return(edge)
}


#' Collapse duplicate INDRA statements into a mapping of edge to metadata
#' @param res INDRA response
#' @param input filtered groupComparison result
#' @importFrom r2r hashmap keys
#' @return processed edge to metadata mapping
#' @keywords internal
#' @noRd
.collapseDuplicateEdgesIntoEdgeToMetadataMapping <- function(res, input) {
    edgeToMetadataMapping <- hashmap()

    for (edge in res) {
        key <- paste(edge$source_id, edge$target_id, edge$data$stmt_type, sep = "_")
        if (key %in% keys(edgeToMetadataMapping)) {
            edgeToMetadataMapping[[key]]$data$evidence_count <-
                edgeToMetadataMapping[[key]]$data$evidence_count +
                edge$data$evidence_count
            edgeToMetadataMapping[[key]]$data$paper_count <- 
                edgeToMetadataMapping[[key]]$data$paper_count + 1
        } else {
            edge <- .addAdditionalMetadataToIndraEdge(edge, input)
            edge$data$paper_count <- 1
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
    # Get unique nodes from edges
    node_ids <- unique(c(edges$source, edges$target))
    
    # Create base nodes dataframe
    nodes <- data.frame(
        id = node_ids,
        stringsAsFactors = FALSE
    )
    
    # Add attributes from input where available
    nodes$logFC <- input$log2FC[match(nodes$id, input$Protein)]
    nodes$adj.pvalue <- input$adj.pvalue[match(nodes$id, input$Protein)]
    nodes$hgncName <- if ("HgncName" %in% colnames(input) && is.character(input$HgncName)) {
        hgnc_value <- input$HgncName[match(nodes$id, input$Protein)]
        ifelse(is.na(hgnc_value), nodes$id, hgnc_value)
    } else {
        nodes$id
    }
    
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
