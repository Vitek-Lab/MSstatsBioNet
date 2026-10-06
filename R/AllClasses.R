#' Network backend classes
#'
#' A backend is a source of prior-knowledge networks, such as INDRA.
#' \code{\link{get_network}()} dispatches on the backend and the query, so
#' each (backend, query) pair has its own method. \code{NetworkBackend} is
#' virtual: create a backend with a constructor such as
#' \code{\link{indra_backend}()}. Other packages can add a backend by
#' extending \code{NetworkBackend} and writing methods for the generics.
#'
#' \code{IndraBackend} queries INDRA CoGEx for networks and grounds names
#' with Gilda, INDRA's grounding service.
#'
#' @slot cogex_url base URL of INDRA CoGEx
#' @slot grounding_url base URL of Gilda
#'
#' @seealso \code{\link{indra_backend}()}, \code{\link{backend_capabilities}()}
#' @name NetworkBackend-class
#' @aliases NetworkBackend-class IndraBackend-class
#' @importFrom methods setClass setValidity
#' @exportClass NetworkBackend IndraBackend
NULL

setClass("NetworkBackend", representation("VIRTUAL"))

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
#' A query says which question \code{\link{get_network}()} asks of the
#' backend. \code{NetworkQuery} is virtual: create a query with a
#' constructor such as \code{\link{subnetwork_query}()}. See
#' \code{\link{network_queries}} for the questions each query answers.
#'
#' @seealso \code{\link{network_queries}}
#' @name NetworkQuery-class
#' @aliases NetworkQuery-class SubnetworkQuery-class
#' @exportClass NetworkQuery SubnetworkQuery
NULL

setClass("NetworkQuery", representation("VIRTUAL"))

setClass("SubnetworkQuery", contains = "NetworkQuery")
