#' Validate input for MSstatsBioNet getSubnetworkFromIndra
#' @param input dataframe from MSstats groupComparison output
#' @param protein_level_data dataframe from MSstats dataProcess output
#' @param sources_filter sources filter
#' @param force_include_other character vector of identifiers to include in the
#' network.
#' @keywords internal
#' @noRd
.validateGetSubnetworkFromIndraInput <- function(input, protein_level_data = NULL, sources_filter, force_include_other) {
    .validateIndraSubnetworkInput(input, sources_filter, force_include_other)
    if (!is.null(protein_level_data)) {
        if(!all(c("Protein", "LogIntensities", "originalRUN") %in% colnames(protein_level_data))) {
            stop("protein_level_data must contain 'Protein', 'LogIntensities', and 'originalRUN' columns.")
        }
    }
}

#' Filter groupComparison result input based on user-defined cutoffs
#' @param input groupComparison result
#' @param pvalueCutoff p-value cutoff
#' @param logfc_cutoff logFC cutoff
#' @param force_include_other list of identifiers to exempt from filtering
#' @param include_infinite_fc logical, whether to include proteins with 
#' infinite log fold change (i.e. proteins that are only detected in one condition).  
#' Default is FALSE. 
#' @param direction Character string specifying the direction of regulation to
#' include. One of \code{"both"} (default), \code{"up"} (upregulated only),
#' or \code{"down"} (downregulated only).
#' @return filtered groupComparison result
#' @keywords internal
#' @noRd
.filterGetSubnetworkFromIndraInput <- function(input,
                                               pvalueCutoff,
                                               logfc_cutoff,
                                               force_include_other,
                                               include_infinite_fc,
                                               direction) {
    input$Protein <- as.character(input$Protein)

    # Drop rows with no grounding (NA EntityId) - they cannot become INDRA
    # query nodes, and keeping them just wastes a fan-out slot.
    if ("EntityId" %in% colnames(input)) {
        na_entity_mask <- is.na(input$EntityId)
        n_dropped <- sum(na_entity_mask)
        if (n_dropped > 0) {
            message("Dropping ", n_dropped,
                    " row(s) with no entity grounding (NA EntityId).")
            input <- input[!na_entity_mask, , drop = FALSE]
        }
    }

    # Extract exempt proteins before any filtering
    exempt_proteins <- NULL
    if (!is.null(force_include_other)) {
        if (!is.character(force_include_other)) {
            stop("force_include_other must be a character vector")
        }
        if ("EntityId" %in% colnames(input) && "EntityNamespace" %in% colnames(input)) {
            fio_pairs <- lapply(force_include_other, function(x) {
                parts <- unlist(strsplit(x, ":"))
                if (length(parts) == 2) list(ns = parts[1], id = parts[2]) else NULL
            })
            fio_pairs <- Filter(Negate(is.null), fio_pairs)
            if (length(fio_pairs) > 0) {
                row_matches_fio <- vapply(seq_len(nrow(input)), function(i) {
                    row_ns <- unlist(strsplit(as.character(input$EntityNamespace[i]), ";"))
                    row_id <- unlist(strsplit(as.character(input$EntityId[i]),        ";"))
                    if (length(row_ns) != length(row_id) || length(row_ns) == 0) {
                        return(FALSE)
                    }
                    for (p in fio_pairs) {
                        if (any(row_ns == p$ns & row_id == p$id)) return(TRUE)
                    }
                    FALSE
                }, logical(1))
                exempt_proteins <- input[row_matches_fio, ]
            } else {
                exempt_proteins <- data.frame()
            }
        } else {
            exempt_proteins <- data.frame()
        }
    }
    
    infinite_fc_proteins <- NULL
    if (include_infinite_fc) {
        infinite_fc_proteins <- input[is.infinite(input$log2FC), ]
    } else {
        if ("log2FC" %in% colnames(input)) {
            input <- input[!is.infinite(input$log2FC), ]
        }
    }

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
    
    if (!is.null(infinite_fc_proteins) && nrow(infinite_fc_proteins) > 0) {
        combined_input <- rbind(infinite_fc_proteins, input)
        input <- combined_input[!duplicated(combined_input$Protein), ]
    }
    
    if (direction == "up") {
        input <- input[!is.na(input$log2FC) & input$log2FC > 0, ]
    } else if (direction == "down") {
        input <- input[!is.na(input$log2FC) & input$log2FC < 0, ]
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

#' Construct nodes data.frame from groupComparison output
#' @param input filtered groupComparison result
#' @param edges edges data frame
#' @return nodes data.frame
#' @keywords internal
#' @noRd
.constructNodesDataFrame <- function(input, edges) {
    nodes = input[, c("Protein", "EntityName", "EntityNamespace", "EntityId",
                      "Site", "log2FC", "adj.pvalue")]
    colnames(nodes) = c("id", "entity_name", "namespace", "entity_id",
                        "site", "log2FC", "adj.pvalue")
    # fread reads numeric IDs (e.g. HGNC) as integers
    for (col in c("id", "entity_name", "namespace", "entity_id", "site")) {
        nodes[[col]] = as.character(nodes[[col]])
    }

    nodes = nodes[nodes$id %in% c(edges$source, edges$target), ]
    extra_force_include_other <- setdiff(unique(c(edges$source, edges$target)), nodes$id)
    if (length(extra_force_include_other) > 0) {
        extra_nodes <- data.frame(
            id = extra_force_include_other,
            entity_name = NA_character_,
            namespace = NA_character_,
            entity_id = NA_character_,
            site = NA_character_,
            log2FC = 0,
            adj.pvalue = 1,
            stringsAsFactors = FALSE
        )
        nodes <- rbind(nodes, extra_nodes)
    }
    nodes$entity_name = ifelse(is.na(nodes$entity_name), nodes$id, nodes$entity_name)

    return(nodes)
}

#' Filter Edges Data Frame
#' @param edges response from INDRA
#' @param correlation_cutoff if protein_level_abundance is not NULL, apply a 
#' cutoff for edges with correlation less than a specified cutoff.
#' Unused when edges have no correlation column.
#' @return filtered edges data frame
#' @keywords internal
#' @noRd
.filterEdgesDataFrame <- function(edges, correlation_cutoff = NULL) {
    if ("correlation" %in% colnames(edges)) {
        edges <- edges[which(abs(edges$correlation) >= correlation_cutoff), ]
    }
    if (nrow(edges) == 0) {
        stop("No edges remain after applying filters. Consider relaxing filters")
    }
    return(edges)
}

#' @importFrom httr GET status_code content
#' @importFrom jsonlite fromJSON
.get_incorrect_curation_count <- function(stmt_hash) {
    stmt_hash_char <- as.character(stmt_hash)
    url <- paste0("https://db.indra.bio/curation/list/", stmt_hash_char)

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

#' Subtract evidence curated as incorrect in INDRA, then re-apply the cutoff
#' @param nodes nodes data frame
#' @param edges edges data frame
#' @param evidence_count_cutoff minimum evidence count per edge
#' @param filter_by_curation logical; if FALSE, nodes and edges are returned
#' unchanged
#' @return list of nodes and edges
#' @keywords internal
#' @noRd
.filterByCuration = function(nodes, edges, evidence_count_cutoff, filter_by_curation) {
    if (filter_by_curation) {
        incorrect_counts <- numeric(nrow(edges))
        for (i in seq_len(nrow(edges))) {
            incorrect_counts[i] <- .get_incorrect_curation_count(edges$statement_id[i])
            Sys.sleep(0.1)
        }
        edges$evidence_count <- as.integer(edges$evidence_count - incorrect_counts)
        edges <- edges[edges$evidence_count >= evidence_count_cutoff, ]
        nodes <- nodes[nodes$id %in% c(edges$source, edges$target), ]
    }
    return(list(nodes = nodes, edges = edges))
}

.filterByPtmSite = function(nodes, edges, filter_by_ptm_site) {
    if (filter_by_ptm_site && nrow(nodes[!is.na(nodes$site), ]) > 0) {
        ptm_overlap <- .ptmOverlap(edges, nodes)
        keep <- ptm_overlap[paste(edges$source, edges$target, edges$interaction, sep = "-")]
        edges <- edges[!is.na(keep) & keep != "", ]
        edges <- edges[!is.na(edges$site),]
        nodes <- nodes[nodes$id %in% c(edges$source, edges$target), ]
    }
    return(list(nodes = nodes, edges = edges))
}

#' Add the correlation between each edge's endpoints
#' @param edges edges data frame
#' @param protein_level_data output of dataProcess
#' @return edges with a correlation column
#' @keywords internal
#' @noRd
.addCorrelationToEdges <- function(edges, protein_level_data) {
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


