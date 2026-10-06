#' Questions you can ask of a network backend
#'
#' Each query constructor asks one question of a backend. Pass the query to
#' \code{\link{get_network}()}. Every query returns the same \code{nodes}
#' and \code{edges} tables (see \code{\link{validate_network}()}), so
#' visualization and filtering work the same whichever question produced
#' the network. The \code{query_type} column of \code{edges} names the
#' query, without \code{_query}.
#'
#' \tabular{llll}{
#'   \strong{Constructor} \tab \strong{Question} \tab
#'   \strong{Nodes added} \tab \strong{\code{query_type}} \cr
#'   \code{\link{subnetwork_query}()} \tab How are my selected entities
#'   connected to each other, with no other nodes added? \tab none, other
#'   than \code{include_entities} \tab \code{"subnetwork"} \cr
#' }
#'
#' More queries (shared regulators, regulator enrichment, paths, and
#' others) are planned. \code{backend_capabilities(backend)$query_types} lists the
#' ones a backend supports.
#'
#' @section Glossary:
#' \describe{
#'   \item{entity}{One row of the entity table from
#'     \code{\link{prepare_entities}()}: one analyte of the input, such as a
#'     protein, a PTM site, or a metabolite.}
#'   \item{selected}{An entity with \code{included_in_query = TRUE}, set by
#'     \code{\link{select_entities}()}: it passed the cutoffs, or was forced
#'     in. Only selected entities are sent to the backend.}
#'   \item{in the input (\code{measured})}{A node that matches a row of the
#'     entity table, whether or not it was selected. It carries that row's
#'     statistics.}
#'   \item{latent}{A node that matches no row of the entity table, such as
#'     an entity added through \code{include_entities}. It has
#'     \code{measured = FALSE} and \code{NA} statistics. It may still have
#'     been measured in another experiment, or filtered out before the
#'     entity table was built.}
#'   \item{grounding}{A \code{"namespace:identifier"} pair that names an
#'     entity in a backend, such as \code{"HGNC:11998"} (TP53).}
#' }
#'
#' @seealso \code{\link{get_network}()}, \code{\link{backend_capabilities}()}
#' @name network_queries
#' @aliases network_queries
NULL

#' How are my selected entities connected to each other?
#'
#' Asks which edges connect the selected entities directly, with no other
#' nodes added. Only entities passed to \code{get_network()} as
#' \code{include_entities} are added.
#'
#' Uses the rows of \code{entities} with \code{included_in_query = TRUE}.
#' Each node the backend returns is matched against all rows, so that it
#' carries its statistics. Edges get \code{query_type = "subnetwork"}.
#' For INDRA, the query goes to CoGEx \code{indra_subnetwork_relations}, and
#' takes fewer than 400 groundings
#' (\code{backend_capabilities(indra_backend())$max_nodes}).
#'
#' @return a \code{SubnetworkQuery} object, to pass to
#' \code{\link{get_network}()}
#' @seealso \code{\link{network_queries}} for the other questions
#' @importFrom methods new
#' @export
#' @examples
#' subnetwork_query()
subnetwork_query <- function() {
    new("SubnetworkQuery")
}
