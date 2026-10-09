# Network provenance, and saving and loading networks.

#' Columns of network$provenance, in order
#' @keywords internal
#' @noRd
PROVENANCE_COLUMNS <- c("backend_database", "query_type", "retrieved_at",
                        "backend_version", "backend_url", "organism",
                        "parameters", "package_version")

#' The current time, in UTC
#'
#' A function of its own so that tests can fix the time.
#' @return POSIXct
#' @keywords internal
#' @noRd
.get_current_time <- function() {
    as.POSIXct(format(Sys.time(), tz = "UTC"), tz = "UTC")
}

#' Build the provenance row of one get_network() call
#'
#' Records when and from where a network was retrieved. Neither INDRA CoGEx
#' nor a live STRING host can return an earlier answer again, so this is
#' the only record of which data a network came from.
#' @param backend_database the \code{edges$backend_database} value
#' @param query_type the \code{edges$query_type} value
#' @param backend_url the URL the query was sent to
#' @param organism NCBI taxon IDs of the queried entity rows
#' @param parameters named list of the query arguments
#' @param backend_version the backend's data version, \code{NA} when the
#' backend has none
#' @param retrieved_at when the response arrived
#' @return data.frame with one row and the columns in
#' \code{PROVENANCE_COLUMNS}
#' @keywords internal
#' @noRd
.build_provenance <- function(backend_database, query_type, backend_url,
                              organism, parameters,
                              backend_version = NA_character_,
                              retrieved_at = .get_current_time()) {
    organism <- unique(organism[!is.na(organism)])
    data.frame(
        backend_database = backend_database,
        query_type       = query_type,
        retrieved_at     = retrieved_at,
        backend_version  = as.character(backend_version),
        backend_url      = backend_url,
        organism         = if (length(organism) == 0) NA_character_ else
            paste(organism, collapse = ";"),
        parameters       = as.character(jsonlite::toJSON(
            parameters, auto_unbox = TRUE, null = "null", digits = NA)),
        package_version  = as.character(
            utils::packageVersion("MSstatsBioNet")),
        stringsAsFactors = FALSE
    )
}

#' Save a network to a file
#'
#' Saves the whole network, including \code{provenance} (when and from
#' which backend it was retrieved) and \code{regulators}, to an
#' \code{.rds} file that \code{\link{load_network}()} reads back.
#' Backends such as INDRA have no data versions, so the same query can
#' give a different network later: save the network to keep the one an
#' analysis used.
#'
#' The network is checked with \code{\link{validate_network}()} first. A
#' network without \code{provenance} is saved with a warning: it is lost
#' when a network is rebuilt with \code{list(nodes = , edges = )}.
#'
#' @param network list of \code{nodes} and \code{edges}, e.g. from
#' \code{\link{get_network}()} or \code{\link{merge_networks}()}
#' @param file path of the file to write, ending in \code{.rds}
#' @return \code{file}, invisibly
#' @seealso \code{\link{load_network}()}
#' @export
#' @examples
#' network <- list(
#'     nodes = data.frame(
#'         id = c("CHK1_HUMAN", "CDC25A_HUMAN"),
#'         entity_type = "protein",
#'         entity_name = c("CHEK1", "CDC25A"),
#'         namespace = "HGNC",
#'         entity_id = c("1925", "1725"),
#'         measured = TRUE,
#'         included_in_query = TRUE,
#'         node_role = "passed_cutoffs"
#'     ),
#'     edges = data.frame(
#'         source = "CHK1_HUMAN",
#'         target = "CDC25A_HUMAN",
#'         interaction = "Phosphorylation",
#'         directed = TRUE,
#'         site = "S76",
#'         confidence = 0.99,
#'         evidence_count = 12L,
#'         evidence_url = paste0("https://db.indra.bio/statements/",
#'                                 "from_hash/-1234?format=html"),
#'         statement_id = "-1234",
#'         backend_database = "INDRA",
#'         query_type = "subnetwork"
#'     ),
#'     provenance = data.frame(
#'         backend_database = "INDRA",
#'         query_type = "subnetwork",
#'         retrieved_at = as.POSIXct("2026-10-09 12:00:00", tz = "UTC")
#'     )
#' )
#' file <- tempfile(fileext = ".rds")
#' save_network(network, file)
#' identical(load_network(file), network)
save_network <- function(network, file) {
    if (!is.character(file) || length(file) != 1 || is.na(file) ||
        !grepl("\\.rds$", file, ignore.case = TRUE)) {
        stop("`file` must be a single path ending in .rds.", call. = FALSE)
    }
    validate_network(network)
    if (is.null(network$provenance)) {
        warning("The network has no provenance, so the file won't record ",
                "when or from which backend it was retrieved. Was it ",
                "rebuilt with list(nodes = , edges = )? Keep the network ",
                "from get_network() or merge_networks() to keep it.",
                call. = FALSE)
    }
    saveRDS(network, file)
    invisible(file)
}

#' Load a network saved with save_network()
#'
#' Reads the network and checks it with \code{\link{validate_network}()},
#' so a file from an older version of the contract gives a clear error
#' here instead of a failure in a later function.
#'
#' @param file path of an \code{.rds} file written by
#' \code{\link{save_network}()}
#' @return the network: list of \code{nodes} and \code{edges}, plus
#' \code{provenance} and \code{regulators} when it had them
#' @seealso \code{\link{save_network}()}
#' @export
#' @examples
#' network <- list(
#'     nodes = data.frame(
#'         id = c("A", "B"), entity_type = "protein",
#'         entity_name = c("A", "B"), namespace = "HGNC",
#'         entity_id = c("1", "2"), measured = TRUE,
#'         included_in_query = TRUE, node_role = "passed_cutoffs"
#'     ),
#'     edges = data.frame(
#'         source = "A", target = "B", interaction = "Activation",
#'         directed = TRUE, site = NA_character_, confidence = 0.9,
#'         evidence_count = 3L, evidence_url = "https://example.org/1",
#'         statement_id = "1", backend_database = "INDRA",
#'         query_type = "subnetwork"
#'     )
#' )
#' file <- tempfile(fileext = ".rds")
#' suppressWarnings(save_network(network, file))
#' load_network(file)$edges
load_network <- function(file) {
    if (!is.character(file) || length(file) != 1 || is.na(file)) {
        stop("`file` must be a single path.", call. = FALSE)
    }
    if (!file.exists(file)) {
        stop("File not found: ", file, call. = FALSE)
    }
    network <- readRDS(file)
    tryCatch(validate_network(network), error = function(e) {
        stop("The network in ", file, " does not meet the current contract. ",
             conditionMessage(e), call. = FALSE)
    })
    network
}
