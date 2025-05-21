#' Get pathways ranked on relevance from INDRA DB
#'
#' @importFrom httr GET content
#' @importFrom MASS fitdistr
#' @importFrom r2r hashmap keys query
#'
#' @param annotated_df output of \code{\link[MSstats]{groupComparison}} function's 
#' comparisionResult table, which contains a list of proteins and their 
#' corresponding p-values, logFCs, along with additional HGNC ID and HGNC 
#' name columns
#' @param main_target A main target, e.g. main target of a drug or protein of
#' particular interest
#' @param target_type One of either 'Protein' or 'Drug'.  Default is 'Protein'
#'
#' @return df of pathways
#'
#' @export
#'
#' @examples
#' annotated_df <- data.table::fread(system.file(
#'     "extdata/groupComparisonModel.csv",
#'     package = "MSstatsBioNet"
#' ))
#' pathways <- getPathwaysFromIndra(annotated_df, "P05067")
#' head(pathways)
#'
getPathwaysFromIndra <- function(annotated_df, main_target = 'MEN1_HUMAN', target_type = "Protein") {
    annotated_df$Protein <- as.character(annotated_df$Protein)
    log2fc_values <- annotated_df$log2FC
    fit <- fitdistr(log2fc_values, "normal")
    para <- fit$estimate
    # pnorm(0, mean = para[1], sd = para[2])
    
    # Fit a power-law-like curve
    b = -0.08
    m = -1.33
    # 1 - 10^(m*log10(1000)+b)
    # probability
    
    # Call INDRA
    if (target_type == "Protein") {
        main_target_row = annotated_df[annotated_df$Protein == main_target,]
        source_id = as.character(main_target_row$HgncId)
        namespace = "@HGNC"
        id_field = "HGNC"
    } else if (target_type == "Drug") {
        source_id = main_target
        namespace = ""
        id_field = "TEXT"
    } else {
        stop("Invalid target type.")
    }
    url = paste('https://db.indra.bio/statements/from_agents?source_idect=',
                source_id, namespace, sep = "")
    response <- GET(url)
    z = content(response)

    edgeToMetadataMapping <- hashmap()
    
    if (length(z$statements) == 0) {
        return(data.frame())
    }
    
    for (index in seq(1, length(z$statements))) {
        edge <- z$statements[[index]]
        if (edge$type == "Complex") {
            if (length(edge$members) == 2) {
                if (identical(edge$members[[1]]$db_refs[[id_field]], source_id)) {
                    obj = edge$members[[2]]$db_refs$HGNC
                    namespaces = names(edge$members[[2]]$db_refs)
                } else {
                    obj = edge$members[[1]]$db_refs$HGNC
                    namespaces = names(edge$members[[1]]$db_refs)
                }
            } else {
                namespaces = c()
            }
        } else if (edge$type == "Phosphorylation") {
            obj = edge$sub$db_refs$HGNC
            namespaces = names(edge$sub$db_refs)
        } else {
            obj = edge$obj$db_refs$HGNC
            namespaces = names(edge$obj$db_refs)
        }
        
        # Filter out edges with no HGNC ID or not in the dataset
        if (!("HGNC" %in% namespaces)) {
            next
        } else if (!(obj %in% annotated_df$HgncId)) {
            next
        }
        
        key <- paste(source_id, obj, sep = "_")
        
        if (key %in% keys(edgeToMetadataMapping)) {
            edgeToMetadataMapping[[key]]$data$evidence_count <-
                edgeToMetadataMapping[[key]]$data$evidence_count +
                z$evidence_counts[[index]]
            edgeToMetadataMapping[[key]]$data$stmt_type <- unique(c(
                edgeToMetadataMapping[[key]]$data$stmt_type,
                edge$type))
        } else {
            # edge <- MSstatsBioNet:::.addAdditionalMetadataToIndraEdge(edge, annotated_df)
            edgeToMetadataMapping[[key]] <- edge
            edgeToMetadataMapping[[key]]$data$evidence_count <-
                z$evidence_counts[[index]]
            edgeToMetadataMapping[[key]]$data$stmt_type <- c(edge$type)
            edgeToMetadataMapping[[key]]$source_id <- source_id
            edgeToMetadataMapping[[key]]$target_id <- obj
            edgeToMetadataMapping[[key]] <- MSstatsBioNet:::.addAdditionalMetadataToIndraEdge(
                edgeToMetadataMapping[[key]], annotated_df, namespace
            )
        }
    }
    
    # Calculate probabilities
    for (key in keys(edgeToMetadataMapping)) {
        edgeToMetadataMapping[[key]]$data$stmt_type <-
            paste(unique(edgeToMetadataMapping[[key]]$data$stmt_type), 
                  collapse = ", ")
        prob_logFC = 0
        logFC = annotated_df[which(annotated_df$HgncId == edgeToMetadataMapping[[key]]$target_id),]
        logFC = logFC$log2FC[[1]]
        if (logFC > para[1]) {
            prob_logFC = 1 - pnorm(logFC, mean = para[1], sd = para[2])
        } else {
            prob_logFC = pnorm(logFC, mean = para[1], sd = para[2])
        }
        evidence_prob = 10^(m*log10(edgeToMetadataMapping[[key]]$data$evidence_count)+b)
        edgeToMetadataMapping[[key]]$data$total_prob = 1 - ((1 - prob_logFC) * (1 - evidence_prob))
        edgeToMetadataMapping[[key]]$data$logFC = logFC
    }
    
    # Construct DF and sort
    edges <- data.frame(
        source = vapply(keys(edgeToMetadataMapping), function(x) {
            query(edgeToMetadataMapping, x)$source_uniprot_id
        }, ""),
        target = vapply(keys(edgeToMetadataMapping), function(x) {
            query(edgeToMetadataMapping, x)$target_uniprot_id
        }, ""),
        interaction = vapply(keys(edgeToMetadataMapping), function(x) {
            query(edgeToMetadataMapping, x)$data$stmt_type
        }, ""),
        evidenceCount = vapply(keys(edgeToMetadataMapping), function(x) {
            query(edgeToMetadataMapping, x)$data$evidence_count
        }, 1),
        logFC = vapply(keys(edgeToMetadataMapping), function(x) {
            query(edgeToMetadataMapping, x)$data$logFC
        }, 1),
        prob = vapply(keys(edgeToMetadataMapping), function(x) {
            query(edgeToMetadataMapping, x)$data$total_prob
        }, 1),
        evidenceLink = vapply(keys(edgeToMetadataMapping), function(x) {
            query(edgeToMetadataMapping, x)$evidence_list
        }, ""),
        stringsAsFactors = FALSE
    )
    
    nodes <- .constructNodesDataFrame(
        annotated_df, edges
    )
    if (!(main_target %in% nodes$id)) {
        # add a row with the main target
        nodes <- rbind(
            nodes,
            data.frame(
                id = main_target,
                logFC = 0,
                pvalue = 0,
                hgncName = main_target
            )
        )
    }
    return(list(nodes = nodes, edges = edges))
}