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
#' statement type in INDRA.  Default is NULL.
#' @param paper_count_cutoff number of papers to filter on. Default is 1.
#' @param evidence_count_cutoff number of evidence to filter on for each
#' paper. E.g. A paper may have 5 sentences describing the same interaction vs 1
#' sentence.  Default is 1.
#' @param correlation_cutoff if protein_level_abundance is not NULL, apply a 
#' cutoff for edges with correlation less than a specified cutoff.  Default is
#' 0.3
#' @param sources_filter filtering only on specific sources.  Default is no filter, i.e. NULL.
#' Otherwise, should be a list, e.g. c('reach', 'medscan').
#' @param logfc_cutoff absolute log fold change cutoff for filtering proteins. 
#' Only proteins with |logFC| greater than this value will be retained. Default 
#' is NULL, i.e. no logFC filtering.
#' @param force_include_other character vector of identifiers to include in the
#' network, regardless if those ids are in the input data. Should be formatted
#' as "namespace:identifier", e.g. "HGNC:1234" or "CHEBI:4911".
#' @param filter_by_curation logical, whether to filter out statements that
#' have been curated as incorrect in INDRA.  Default is FALSE.
#' @param filter_by_ptm_site logical, whether to filter edges based on whether the 
#' site information from INDRA matches with the PTM site in the input.  Default is FALSE.  
#' Only applicable for differential PTM abundance results.
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
                                   statement_types = NULL,
                                   paper_count_cutoff = 1,
                                   evidence_count_cutoff = 1,
                                   correlation_cutoff = 0.3,
                                   sources_filter = NULL,
                                   logfc_cutoff = NULL,
                                   force_include_other = NULL, 
                                   filter_by_curation = FALSE,
                                   filter_by_ptm_site = FALSE) {
    input <- .filterGetSubnetworkFromIndraInput(input, pvalueCutoff, logfc_cutoff, force_include_other)
    .validateGetSubnetworkFromIndraInput(input, protein_level_data, sources_filter, force_include_other)
    res <- .callIndraCogexApi(input$HgncId, force_include_other)
    res <- .filterIndraResponse(res, statement_types, evidence_count_cutoff, sources_filter)
    edges <- .constructEdgesDataFrame(res, input, protein_level_data)
    edges <- .filterEdgesDataFrame(edges, paper_count_cutoff, correlation_cutoff)
    nodes <- .constructNodesDataFrame(input, edges)
    subnetwork = .filterByPtmSite(nodes, edges, filter_by_ptm_site)
    subnetwork = .filterByCuration(subnetwork$nodes, subnetwork$edges, filter_by_curation)
    warning(
        "NOTICE: This function includes third-party software components
        that are licensed under the BSD 2-Clause License. Please ensure to
        include the third-party licensing agreements if redistributing this
        package or utilizing the results based on this package.
        See the LICENSE file for more details."
    )
    return(subnetwork)
}
