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
getPathwaysFromIndra <- function(annotated_df, main_target = 'MEN1_HUMAN') {
    log2fc_values <- annotated_df$log2FC
    fit <- fitdistr(log2fc_values, "normal")
    para <- fit$estimate
    # pnorm(0, mean = para[1], sd = para[2])
    
    # Fit a negative binomial distribution with parameters
    n = 0.5
    p = 0.06
    # probability <- dnbinom(2, size = n, prob = p)
    # probability
    
    # Call INDRA
    main_target_row = annotated_df[annotated_df$Protein == main_target,]
    source_id = main_target_row$HgncId
    url = paste('https://db.indra.bio/statements/from_agents?subject=',
                source_id, '@HGNC', sep = "")
    response <- GET(url)
    z = content(response)

    edgeToMetadataMapping <- hashmap()
    
    if (length(z$statements) == 0) {
        return(data.frame())
    }
    
    for (index in seq(1, length(z$statements))) {
        edge <- z$statements[[index]]
        key <- ""
        if (edge$type == "Complex") {
            next
        } else if (!("HGNC" %in% names(edge$obj$db_refs))) {
            next
        } else if (!(edge$obj$db_refs$HGNC %in% annotated_df$HgncId)) {
            next
        } else {
            key <- paste(edge$subj$db_refs$HGNC, edge$obj$db_refs$HGNC, sep = "_")
        }
        
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
            edgeToMetadataMapping[[key]]$source_id <- edge$subj$db_refs$HGNC
            edgeToMetadataMapping[[key]]$target_id <- edge$obj$db_refs$HGNC
            edgeToMetadataMapping[[key]] <- MSstatsBioNet:::.addAdditionalMetadataToIndraEdge(
                edgeToMetadataMapping[[key]], annotated_df
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
        evidence_prob = dnbinom(min(10, edgeToMetadataMapping[[key]]$data$evidence_count), size = n, prob = p)
        edgeToMetadataMapping[[key]]$data$total_prob = prob_logFC * evidence_prob
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
}