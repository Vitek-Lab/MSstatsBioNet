# Filters that run on a finished network, after get_network().

#' Remove evidence curated as incorrect
#'
#' Subtracts the evidence curated as incorrect, from
#' \code{\link{get_curations}()}, from each edge's \code{evidence_count},
#' then drops the edges left below \code{min_evidence}, and the nodes those
#' edges leave without any edge. Edges left with no evidence are dropped
#' whatever \code{min_evidence} is. A message says how many edges were
#' dropped.
#'
#' Run it after any filter that pools edges, such as the PTM-site filter of
#' \code{\link{getSubnetworkFromIndra}()}: dropping edges first can change
#' which edges such a filter keeps.
#'
#' @param network list of \code{nodes} and \code{edges} from
#' \code{\link{get_network}()}. \code{edges} needs the columns
#' \code{source}, \code{target}, \code{statement_id}, and
#' \code{evidence_count}, and \code{backend_database} unless
#' \code{backend} is given.
#' @param min_evidence minimum evidence count per edge, after subtracting
#' the incorrect evidence
#' @param backend the backend to get the curations from, e.g.
#' \code{\link{indra_backend}()}. \code{NULL} (default) uses the default
#' backend named in each edge's \code{backend_database}. Pass a backend to
#' use one built with non-default settings, such as
#' \code{indra_backend(curation_url = )}. A backend without curations, such
#' as one with no \code{get_curations()} method, gives an error naming it.
#' @return \code{network}, with \code{edges$evidence_count} reduced and the
#' dropped edges and nodes removed. Other elements of \code{network} are
#' kept.
#' @seealso \code{\link{get_curations}()}
#' @export
#' @examples
#' \donttest{
#' input <- data.table::fread(system.file(
#'     "extdata/groupComparisonModel.csv",
#'     package = "MSstatsBioNet"
#' ))
#' indra <- indra_backend()
#' entities <- prepare_entities(input, entity_type = "protein",
#'                              id_type = "uniprot")
#' entities <- convert_ids(indra, entities)
#' entities <- select_entities(entities, pvalue_cutoff = 0.05)
#' network <- get_network(indra, entities, interaction_types = "Complex")
#' network <- filter_by_curation(network, min_evidence = 2)
#' head(network$edges)
#' }
filter_by_curation <- function(network, min_evidence = 1, backend = NULL) {
    .check_filter_by_curation_input(network, min_evidence)
    .check_backend_argument(backend)
    edges <- network$edges
    incorrect_counts <- .count_incorrect_evidence(edges, backend)
    edges$evidence_count <- as.integer(edges$evidence_count - incorrect_counts)
    keep <- edges$evidence_count >= min_evidence & edges$evidence_count >= 1
    if (any(!keep)) {
        message("Dropping ", sum(!keep), " edge(s) with fewer than ",
                max(min_evidence, 1), " evidence after removing the ",
                "evidence curated as incorrect.")
    }
    nodes <- network$nodes
    endpoints_before <- c(edges$source, edges$target)
    edges <- edges[keep, , drop = FALSE]
    endpoints_after <- c(edges$source, edges$target)
    left_without_edges <- nodes$id %in% endpoints_before &
        !nodes$id %in% endpoints_after
    network$nodes <- nodes[!left_without_edges, , drop = FALSE]
    network$edges <- edges
    network
}

#' Check the input of filter_by_curation()
#' @param network the network argument
#' @param min_evidence the min_evidence argument
#' @return \code{NULL}, invisibly; errors otherwise
#' @keywords internal
#' @noRd
.check_filter_by_curation_input <- function(network, min_evidence) {
    if (!is.list(network) || !is.data.frame(network$nodes) ||
        !is.data.frame(network$edges)) {
        stop("`network` must be a list with `nodes` and `edges` ",
             "data.frames, e.g. from get_network().", call. = FALSE)
    }
    missing_cols <- setdiff(c("source", "target", "statement_id",
                              "evidence_count"), names(network$edges))
    if (length(missing_cols) > 0) {
        stop("`network$edges` is missing the column(s): ",
             paste(missing_cols, collapse = ", "), call. = FALSE)
    }
    if (!"id" %in% names(network$nodes)) {
        stop("`network$nodes` is missing the column: id", call. = FALSE)
    }
    if (!is.numeric(min_evidence) || length(min_evidence) != 1 ||
        !is.finite(min_evidence)) {
        stop("`min_evidence` must be a single number.", call. = FALSE)
    }
    invisible(NULL)
}
