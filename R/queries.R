#' Subnetwork query
#'
#' Asks which edges connect the selected entities directly. It adds no
#' nodes, other than entities passed as \code{include_entities} to
#' \code{get_network()}. Edges get \code{query_type = "subnetwork"}.
#'
#' Internal until the entity model is added (Phase 3 of the API refactor).
#' @return a \code{SubnetworkQuery} object
#' @importFrom methods new
#' @keywords internal
#' @noRd
subnetwork_query <- function() {
    new("SubnetworkQuery")
}
