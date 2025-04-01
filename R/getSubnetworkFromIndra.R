#' Get subnetwork from INDRA database
#'
#' Using differential abundance results from MSstats, this function retrieves
#' a subnetwork of protein interactions from INDRA database.
#'
#' @param input output of \code{\link[MSstats]{groupComparison}} function's 
#' comparisionResult table, which contains a list of proteins and their 
#' corresponding p-values, logFCs, along with additional HGNC ID and HGNC 
#' name columns
#' @param protein_level_data output of the \code{\link[MSstats]{dataProcess}} 
#' function's ProteinLevelData table, which contains a list of proteins and 
#' their corresponding abundances.  Used for annotating correlation information 
#' and applying correlation cutoffs.
#' @param pvalueCutoff p-value cutoff for filtering. Default is NULL, i.e. no
#' filtering
#' @param statement_types list of interaction types to filter on.  Equivalent to
#' statement type in INDRA.  Default is c("IncreaseAmount", "DecreaseAmount").
#' @param paper_count_cutoff number of papers to filter on. Default is 1.
#' @param evidence_count_cutoff number of evidence to filter on for each
#' paper. E.g. A paper may have 5 sentences describing the same interaction vs 1
#' sentence.  Default is 1.
#' @param correlation_cutoff if protein_level_abundance is not NULL, apply a 
#' cutoff for edges with correlation less than a specified cutoff.  Default is
#' 0.3
#' @param sources_filter filtering only on specific sources.  Default is no filter, i.e. NULL.
#' Otherwise, should be a list, e.g. c('reach', 'medscan').
#'
#' @return list of 2 data.frames, nodes and edges
#'
#' @export
#'
#' @examples
#' input <- data.table::fread(system.file(
#'     "extdata/groupComparisonModel.csv",
#'     package = "MSstatsBioNet"
#' ))
#' subnetwork <- getSubnetworkFromIndra(input)
#' head(subnetwork$nodes)
#' head(subnetwork$edges)
#'
getSubnetworkFromIndra <- function(input, 
                                   protein_level_data = NULL,
                                   pvalueCutoff = NULL, 
                                   statement_types = c("IncreaseAmount", "DecreaseAmount"),
                                   paper_count_cutoff = 1,
                                   evidence_count_cutoff = 1,
                                   correlation_cutoff = 0.3,
                                   sources_filter = NULL) {
    input <- .filterGetSubnetworkFromIndraInput(input, pvalueCutoff)
    .validateGetSubnetworkFromIndraInput(input, protein_level_data)
    res <- .callIndraCogexApi(input$HgncId)
    res <- .filterIndraResponse(res, statement_types, evidence_count_cutoff, sources_filter)
    edges <- .constructEdgesDataFrame(res, input, protein_level_data)
    edges <- .filterEdgesDataFrame(edges, paper_count_cutoff, correlation_cutoff)
    nodes <- .constructNodesDataFrame(input, edges)
    warning(
        "NOTICE: This function includes third-party software components
        that are licensed under the BSD 2-Clause License. Please ensure to
        include the third-party licensing agreements if redistributing this
        package or utilizing the results based on this package.
        See the LICENSE file for more details."
    )
    return(list(nodes = nodes, edges = edges))
}
