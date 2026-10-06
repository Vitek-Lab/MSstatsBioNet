# The input filter getSubnetworkFromIndra() used before Phase 3c of the API
# refactor, kept as the reference for the select_entities() equivalence test
# in test-entities.R. It drops the rows that fail the cutoffs, where
# select_entities() flags them.
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
