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
#' with Gilda, INDRA's grounding service. Curations come from the INDRA
#' database.
#'
#' \code{StringBackend} queries the STRING database of protein-protein
#' associations through its REST API, for one network type and one STRING
#' version.
#'
#' @slot cogex_url (\code{IndraBackend}) base URL of INDRA CoGEx
#' @slot grounding_url (\code{IndraBackend}) base URL of Gilda
#' @slot curation_url (\code{IndraBackend}) base URL of the INDRA
#'   database, which holds the curations
#' @slot network_type (\code{StringBackend}) \code{"physical"},
#'   \code{"functional"}, or \code{"regulatory"}
#' @slot version (\code{StringBackend}) the STRING version, e.g.
#'   \code{"12.5"}, or \code{NA} for the current one
#' @slot caller_identity (\code{StringBackend}) the name sent to STRING
#'   with each request
#' @slot host (\code{StringBackend}) environment holding the STRING
#'   version and address, once resolved
#'
#' @seealso \code{\link{indra_backend}()}, \code{\link{string_backend}()},
#' \code{\link{backend_capabilities}()}
#' @name NetworkBackend-class
#' @aliases NetworkBackend-class IndraBackend-class StringBackend-class
#' @importFrom methods setClass setValidity
#' @exportClass NetworkBackend IndraBackend StringBackend
NULL

setClass("NetworkBackend", representation("VIRTUAL"))

setClass("IndraBackend", contains = "NetworkBackend",
         representation(cogex_url = "character",
                        grounding_url = "character",
                        curation_url = "character"))

setValidity("IndraBackend", function(object) {
    for (slot_name in c("cogex_url", "grounding_url", "curation_url")) {
        url <- methods::slot(object, slot_name)
        if (length(url) != 1 || is.na(url) || !nzchar(url)) {
            return(paste(slot_name, "must be a single non-empty string"))
        }
    }
    TRUE
})

setClass("StringBackend", contains = "NetworkBackend",
         representation(network_type = "character",
                        version = "character",
                        caller_identity = "character",
                        host = "environment"))

setValidity("StringBackend", function(object) {
    if (length(object@network_type) != 1 ||
        !object@network_type %in% STRING_NETWORK_TYPES) {
        return(paste0("network_type must be one of ",
                      paste0("\"", STRING_NETWORK_TYPES, "\"",
                             collapse = ", ")))
    }
    if (length(object@version) != 1 ||
        (!is.na(object@version) &&
         !grepl("^[0-9]+\\.[0-9]+$", object@version))) {
        return("version must be NULL or a single string such as \"12.5\"")
    }
    if (length(object@caller_identity) != 1 ||
        is.na(object@caller_identity) || !nzchar(object@caller_identity)) {
        return("caller_identity must be a single non-empty string")
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
