#' Create an INDRA backend
#'
#' Internal until the entity model is added (Phase 3 of the API refactor).
#' @param cogex_url base URL of INDRA CoGEx
#' @return an \code{IndraBackend} object
#' @importFrom methods new
#' @keywords internal
#' @noRd
indra_backend <- function(cogex_url = INDRA_API_URL) {
    new("IndraBackend", cogex_url = cogex_url)
}

#' INDRA subnetwork query
#'
#' Sends the selected rows' groundings to CoGEx
#' \code{indra_subnetwork_relations}, filters the statements, and normalizes
#' them to the edge contract. Nodes are the selected rows that an edge
#' reaches, plus any \code{include_entities} endpoints.
#' @keywords internal
#' @noRd
setMethod("get_network", signature("IndraBackend", "SubnetworkQuery"),
    function(backend, entities, query, statement_types = NULL,
             min_evidence = 1, sources = NULL, include_entities = NULL,
             ...) {
        .validateIndraSubnetworkInput(entities, sources, include_entities)
        res <- .callIndraCogexApi(entities$EntityNamespace, entities$EntityId,
                                  include_entities, backend@cogex_url)
        res <- .filterIndraResponse(res, statement_types, min_evidence,
                                    sources)
        edges <- .constructEdgesDataFrame(res, entities)
        edges <- .filterEdgesDataFrame(edges)
        network <- list(nodes = .constructNodesDataFrame(entities, edges),
                        edges = edges)
        validate_network(network)
        network
    })

#' Validate the input of the INDRA subnetwork query
#' @param input annotated groupComparison table of the selected rows
#' @param sources_filter sources filter
#' @param force_include_other character vector of identifiers to include in
#' the network
#' @keywords internal
#' @noRd
.validateIndraSubnetworkInput <- function(input, sources_filter, force_include_other) {
    required_cols <- c("Protein", "log2FC", "adj.pvalue",
                       "EntityNamespace", "EntityId", "EntityName")
    missing_cols <- setdiff(required_cols, colnames(input))
    if (length(missing_cols) > 0) {
        stop("Invalid Input Error: input is missing required column(s): ",
             paste(missing_cols, collapse = ", "), ".")
    }
    ids_split <- unlist(strsplit(as.character(input$EntityId),        ";"), use.names = FALSE)
    nss_split <- unlist(strsplit(as.character(input$EntityNamespace), ";"), use.names = FALSE)
    unique_pairs <- unique(paste(nss_split, ids_split, sep = ":"))
    num_proteins = length(unique_pairs) +
        ifelse(!is.null(force_include_other), length(force_include_other), 0)
    if (num_proteins >= 400) {
        stop("Invalid Input Error: INDRA query must contain less than 400 proteins.  Consider lowering your p-value cutoff")
    }
    if (nrow(input) == 0) {
        stop("Invalid Input Error: Input must contain at least one protein after filtering.")
    }
    if (!is.null(sources_filter)) {
        if (!is.character(sources_filter)) {
            stop("sources_filter must be a character vector")
        }
    }
}
