#' Get subnetwork from INDRA database
#'
#' Using differential abundance results from MSstats, this function retrieves
#' a subnetwork of protein interactions from INDRA database.
#'
#' @param input output of \code{\link[MSstats]{groupComparison}} function's
#' comparisionResult table, annotated by
#' \code{\link{annotateProteinInfoFromIndra}}. Must contain \code{Protein},
#' \code{EntityNamespace}, and \code{EntityId} columns (and typically also
#' \code{EntityName}, \code{log2FC}, \code{adj.pvalue}). When an analyte
#' grounds to multiple candidates the three \code{Entity*} columns are
#' semicolon-joined and positionally aligned.
#' @param protein_level_data Deprecated, and will be removed in a future
#' release. Output of the \code{\link[MSstats]{dataProcess}} function's
#' ProteinLevelData table, used to annotate edges with correlations and apply
#' \code{correlation_cutoff}. Supplying it gives a deprecation warning.
#' @param pvalueCutoff p-value cutoff for filtering. Default is NULL, i.e. no
#' filtering
#' @param statement_types list of interaction types to filter on.  Equivalent to
#' statement type in INDRA.  Default is NULL.
#' @param paper_count_cutoff Deprecated, and will be removed in a future
#' release. It is ignored: paper counts are not available from INDRA, so this
#' filter never had an effect for 1 and removed every edge for larger values.
#' Supplying it gives a deprecation warning.
#' @param evidence_count_cutoff number of evidence to filter on for each
#' paper. E.g. A paper may have 5 sentences describing the same interaction vs 1
#' sentence.  Default is 1.
#' @param correlation_cutoff Deprecated, and will be removed in a future
#' release. If \code{protein_level_data} is not NULL, remove edges whose
#' absolute correlation is below this cutoff. Default is 0.3. Supplying it
#' gives a deprecation warning.
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
#' @param include_infinite_fc logical, whether to include proteins with 
#' infinite log fold change (i.e. proteins that are only detected in one condition).  
#' Default is FALSE.
#' @param direction Character string specifying the direction of regulation to
#' include. One of \code{"both"} (default), \code{"up"} (upregulated only),
#' or \code{"down"} (downregulated only).
#'
#' @return list of 2 data.frames, \code{nodes} and \code{edges}, that meets
#' the contract checked by \code{\link{validate_network}}.
#'
#' \code{edges} has one row per INDRA statement: \code{source},
#' \code{target}, \code{interaction} (INDRA statement type),
#' \code{directed} (\code{FALSE} for symmetric types such as
#' \code{Complex}), \code{site} (PTM site on the target, or \code{NA}),
#' \code{confidence} (INDRA belief score), \code{evidence_count},
#' \code{evidence_url} (INDRA page for this statement), \code{statement_id}
#' (INDRA statement hash, as character), \code{backend_database}
#' (\code{"INDRA"}), \code{query_type} (\code{"subnetwork"}),
#' \code{evidence_sources} (evidence count per source, as JSON), and the
#' deprecated \code{paperCount} (and \code{correlation} when
#' \code{protein_level_data} is given).
#'
#' \code{nodes} has one row per analyte: \code{id}, \code{entity_name},
#' \code{namespace}, \code{entity_id}, \code{site}, \code{log2FC}, and
#' \code{adj.pvalue}.
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
                                   filter_by_ptm_site = FALSE,
                                   include_infinite_fc = FALSE,
                                   direction = c("both", "up", "down")) {
    direction = match.arg(direction)
    if (!missing(paper_count_cutoff)) {
        .warn_deprecated_arg("paper_count_cutoff", "It is ignored.")
    }
    if (!missing(correlation_cutoff)) {
        .warn_deprecated_arg("correlation_cutoff")
    }
    if (!is.null(protein_level_data)) {
        .warn_deprecated_arg("protein_level_data")
    }
    input <- .filterGetSubnetworkFromIndraInput(input, pvalueCutoff, logfc_cutoff, force_include_other, include_infinite_fc, direction)
    .validateGetSubnetworkFromIndraInput(input, protein_level_data, sources_filter, force_include_other)
    res <- .callIndraCogexApi(input$EntityNamespace, input$EntityId, force_include_other)
    res <- .filterIndraResponse(res, statement_types, evidence_count_cutoff, sources_filter)
    edges <- .constructEdgesDataFrame(res, input, protein_level_data)
    edges <- .filterEdgesDataFrame(edges, correlation_cutoff)
    nodes <- .constructNodesDataFrame(input, edges)
    subnetwork = .filterByPtmSite(nodes, edges, filter_by_ptm_site)
    subnetwork = .filterByCuration(subnetwork$nodes, subnetwork$edges, evidence_count_cutoff, filter_by_curation)
    validate_network(subnetwork)
    warning(
        "NOTICE: This function includes third-party software components
        that are licensed under the BSD 2-Clause License. Please ensure to
        include the third-party licensing agreements if redistributing this
        package or utilizing the results based on this package.
        See the LICENSE file for more details."
    )
    return(subnetwork)
}

#' Warn that an argument of getSubnetworkFromIndra is deprecated
#' @param arg name of the deprecated argument
#' @param detail optional extra sentence describing current behavior
#' @keywords internal
#' @noRd
.warn_deprecated_arg <- function(arg, detail = NULL) {
    warning(
        "Argument '", arg, "' of getSubnetworkFromIndra() is deprecated ",
        "and will be removed in a future release. ",
        if (!is.null(detail)) paste0(detail, " "),
        "See NEWS.md for details.",
        call. = FALSE
    )
}
