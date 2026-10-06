#' Validate input for MSstatsBioNet getSubnetworkFromIndra
#'
#' Checks the input columns. The query size and \code{sources_filter} are
#' checked by \code{get_network()}.
#' @param input dataframe from MSstats groupComparison output
#' @param protein_level_data dataframe from MSstats dataProcess output
#' @param force_include_other character vector of identifiers to include in
#' the network.
#' @keywords internal
#' @noRd
.validateGetSubnetworkFromIndraInput <- function(input, protein_level_data = NULL,
                                                 force_include_other = NULL) {
    required_columns <- c("Protein", "log2FC", "adj.pvalue",
                       "EntityNamespace", "EntityId", "EntityName")
    missing_columns <- setdiff(required_columns, colnames(input))
    if (length(missing_columns) > 0) {
        stop("Invalid Input Error: input is missing required column(s): ",
             paste(missing_columns, collapse = ", "), ".")
    }
    if ("Label" %in% colnames(input)) {
        labels <- unique(input$Label[!is.na(input$Label)])
        if (length(labels) > 1) {
            stop("input has ", length(labels), " comparisons in its Label ",
                 "column: ", .list_values_for_message(labels), ". Filter ",
                 "input to one comparison.", call. = FALSE)
        }
    }
    if (!is.null(force_include_other) && !is.character(force_include_other)) {
        stop("force_include_other must be a character vector")
    }
    if (!is.null(protein_level_data)) {
        if(!all(c("Protein", "LogIntensities", "originalRUN") %in% colnames(protein_level_data))) {
            stop("protein_level_data must contain 'Protein', 'LogIntensities', and 'originalRUN' columns.")
        }
    }
}

#' Build an entity table from input annotated by annotateProteinInfoFromIndra
#'
#' The input is already grounded, so its \code{EntityNamespace},
#' \code{EntityId}, and \code{EntityName} columns are copied to the
#' grounding columns, and \code{convert_ids()} is not needed. Rows whose
#' \code{Protein} ends in a PTM site (\code{_S148}) are \code{ptm_site} rows,
#' with \code{GlobalProtein}, when present, as the parent. Other rows get
#' their entity type from their first namespace (\code{"protein"} when
#' ungrounded). \code{id_type} is not read again, so it is set to the usual
#' identifier system of the entity type.
#' @param input annotated groupComparison result
#' @return entity table
#' @keywords internal
#' @noRd
.build_entities_from_annotated_input <- function(input) {
    input <- as.data.frame(input)
    input$Protein <- as.character(input$Protein)
    first_namespace <- vapply(strsplit(as.character(input$EntityNamespace),
                                       ";", fixed = TRUE),
                              function(namespaces) namespaces[1],
                              character(1))
    entity_types <- .get_namespace_entity_types(first_namespace)
    entity_types[is.na(first_namespace) | entity_types == "other"] <- "protein"
    entity_types[!is.na(parse_ptm_sites(input$Protein)$site)] <- "ptm_site"
    input$inferred_entity_type <- entity_types
    input$inferred_id_type <- ifelse(entity_types %in% c("protein", "ptm_site"),
                             "uniprot", "chemical_name")
    entities <- prepare_entities(input, id_column = "Protein",
                                 entity_type = "inferred_entity_type",
                                 id_type = "inferred_id_type",
                                 logfc_column = "log2FC")
    # fread reads numeric IDs (e.g. HGNC) as integers
    entities$namespace <- as.character(input$EntityNamespace)
    entities$entity_id <- as.character(input$EntityId)
    entities$entity_name <- as.character(input$EntityName)
    entities
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


