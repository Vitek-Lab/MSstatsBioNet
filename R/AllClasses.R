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
#' Queries INDRA CoGEx, and grounds names with Gilda. The Network Search URL
#' is added with its first query (Phase 7 of the API refactor).
#' @slot cogex_url base URL of INDRA CoGEx
#' @slot grounding_url base URL of Gilda
#' @keywords internal
#' @noRd
setClass("IndraBackend", contains = "NetworkBackend",
         representation(cogex_url = "character",
                        grounding_url = "character"))

setValidity("IndraBackend", function(object) {
    for (slot_name in c("cogex_url", "grounding_url")) {
        url <- methods::slot(object, slot_name)
        if (length(url) != 1 || is.na(url) || !nzchar(url)) {
            return(paste(slot_name, "must be a single non-empty string"))
        }
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
