#' Network backend classes
#'
#' A backend is a source of prior-knowledge networks. \code{get_network()}
#' dispatches on the backend and the query, so each (backend, query) pair has
#' its own method.
#'
#' Internal until the entity model is added (Phase 3 of the API refactor).
#' @importFrom methods setClass setValidity
#' @keywords internal
#' @noRd
setClass("NetworkBackend", representation("VIRTUAL"))

#' INDRA backend
#'
#' Queries INDRA CoGEx. The Network Search and Gilda URLs are added with
#' their first queries (Phase 7 and Phase 3 of the API refactor).
#' @slot cogex_url base URL of INDRA CoGEx
#' @keywords internal
#' @noRd
setClass("IndraBackend", contains = "NetworkBackend",
         representation(cogex_url = "character"))

setValidity("IndraBackend", function(object) {
    if (length(object@cogex_url) != 1 || is.na(object@cogex_url) ||
        !nzchar(object@cogex_url)) {
        return("cogex_url must be a single non-empty string")
    }
    TRUE
})

#' Network query classes
#'
#' A query says which question \code{get_network()} asks of the backend.
#' @keywords internal
#' @noRd
setClass("NetworkQuery", representation("VIRTUAL"))

#' Subnetwork query: the edges among the selected entities, adding no nodes
#' @keywords internal
#' @noRd
setClass("SubnetworkQuery", contains = "NetworkQuery")
