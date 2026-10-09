#' Default backend for each backend_database value
#'
#' Used when a function that takes a finished network is called with
#' \code{backend = NULL}: each edge goes to the backend named in its
#' \code{backend_database}, built with default settings.
#' @keywords internal
#' @noRd
BACKEND_CONSTRUCTORS <- list(
    INDRA  = function() indra_backend(),
    STRING = function() string_backend()
)

#' Stop when a backend can't convert some rows of an entity table
#'
#' Checks each row's (\code{entity_type}, \code{id_type}) pair against the
#' backend's supported conversions.
#' @param entities entity table
#' @param backend the backend, for the message
#' @param conversions named list: for each \code{entity_type}, the
#' \code{id_type} values the backend converts
#' @return \code{NULL}, invisibly; errors otherwise
#' @keywords internal
#' @noRd
.check_id_conversions <- function(entities, backend, conversions) {
    pairs <- unique(entities[, c("entity_type", "id_type")])
    supported <- vapply(seq_len(nrow(pairs)), function(i) {
        pairs$id_type[i] %in% conversions[[pairs$entity_type[i]]]
    }, logical(1))
    if (any(!supported)) {
        unsupported <- paste(pairs$entity_type[!supported],
                             pairs$id_type[!supported], sep = " / ")
        allowed <- unlist(lapply(names(conversions), function(type) {
            paste(type, conversions[[type]], sep = " / ")
        }))
        stop(class(backend), " can't convert entity_type / id_type: ",
             .list_values_for_message(unsupported), ". Supported: ",
             paste(allowed, collapse = ", "), ".", call. = FALSE)
    }
    invisible(NULL)
}

#' Check a backend argument
#' @param backend \code{NULL} or a \code{NetworkBackend}
#' @return \code{NULL}, invisibly; errors otherwise
#' @keywords internal
#' @noRd
.check_backend_argument <- function(backend) {
    if (!is.null(backend) && !methods::is(backend, "NetworkBackend")) {
        stop("`backend` must be NULL or a NetworkBackend, e.g. from ",
             "indra_backend().", call. = FALSE)
    }
    invisible(NULL)
}

#' Split edges by the backend that answers for them
#'
#' With \code{backend = NULL}, the edges are split by
#' \code{backend_database} and each group gets that value's default backend
#' from \code{BACKEND_CONSTRUCTORS}. A given \code{backend} takes all edges.
#'
#' @param edges edges data.frame
#' @param backend \code{NULL} or a \code{NetworkBackend}
#' @return list of \code{list(backend =, edges =)}, one per backend, in
#'   order of first appearance in \code{edges}
#' @keywords internal
#' @noRd
.split_edges_by_backend <- function(edges, backend = NULL) {
    .check_backend_argument(backend)
    if (!is.null(backend)) {
        return(list(list(backend = backend, edges = edges)))
    }
    if (!"backend_database" %in% names(edges)) {
        stop("`edges` has no `backend_database` column, so the backend ",
             "can't be chosen. Pass `backend`, e.g. ",
             "`backend = indra_backend()`.", call. = FALSE)
    }
    databases <- as.character(edges$backend_database)
    if (anyNA(databases) || any(!nzchar(databases))) {
        stop("`backend_database` is missing for some edges. Pass ",
             "`backend`, e.g. `backend = indra_backend()`.", call. = FALSE)
    }
    unknown <- setdiff(unique(databases), names(BACKEND_CONSTRUCTORS))
    if (length(unknown) > 0) {
        stop("No default backend for backend_database ",
             paste0("\"", unknown, "\"", collapse = ", "),
             ". Pass `backend` to choose one.", call. = FALSE)
    }
    lapply(unique(databases), function(database) {
        list(backend = BACKEND_CONSTRUCTORS[[database]](),
             edges = edges[databases == database, , drop = FALSE])
    })
}

#' Get the evidence of edges from the backend of each edge
#'
#' Calls \code{get_evidence()} once per backend, see
#' \code{.split_edges_by_backend()}, and stacks the results.
#'
#' @param edges edges data.frame
#' @param backend \code{NULL} or a \code{NetworkBackend}
#' @return evidence data.frame, as from \code{get_evidence()}
#' @keywords internal
#' @noRd
.fetch_evidence <- function(edges, backend = NULL) {
    groups <- .split_edges_by_backend(edges, backend)
    if (length(groups) == 0) {
        warning("No evidence text found for any statement hash")
        return(.build_empty_evidence_table())
    }
    evidence <- lapply(groups, function(group) {
        get_evidence(group$backend, group$edges)
    })
    if (length(evidence) == 1) {
        return(evidence[[1]])
    }
    do.call(rbind, evidence)
}

#' Count the incorrect evidence of each edge, from the backend of each edge
#'
#' Calls \code{get_curations()} once per backend, see
#' \code{.split_edges_by_backend()}, and matches the counts back to the
#' edges within each backend, since two backends can share a
#' \code{statement_id}. Statements a backend returns no count for count 0,
#' after \code{.check_curation_table()} has ruled out a \code{statement_id}
#' that can't match, such as a number.
#'
#' @param edges edges data.frame
#' @param backend \code{NULL} or a \code{NetworkBackend}
#' @return integer vector, one count per row of \code{edges}
#' @keywords internal
#' @noRd
.count_incorrect_evidence <- function(edges, backend = NULL) {
    incorrect_counts <- integer(nrow(edges))
    edges$.edge_row <- seq_len(nrow(edges))
    for (group in .split_edges_by_backend(edges, backend)) {
        curations <- get_curations(group$backend, group$edges)
        .check_curation_table(curations, group$backend)
        counts <- curations$incorrect_count[match(
            as.character(group$edges$statement_id),
            as.character(curations$statement_id))]
        counts[is.na(counts)] <- 0L
        incorrect_counts[group$edges$.edge_row] <- as.integer(counts)
    }
    incorrect_counts
}

#' Check the table a get_curations() method returned
#'
#' A numeric \code{statement_id} matches no edge (INDRA hashes above 2^53
#' lose precision as numbers), so every edge would silently count 0
#' incorrect evidence. Errors naming the backend.
#'
#' @param curations output of \code{get_curations()}
#' @param backend the backend that returned it
#' @return \code{NULL}, invisibly; errors otherwise
#' @keywords internal
#' @noRd
.check_curation_table <- function(curations, backend) {
    problem <- if (!is.data.frame(curations) ||
                   !all(c("statement_id", "incorrect_count") %in%
                        names(curations))) {
        "a data.frame with columns statement_id and incorrect_count"
    } else if (!is.character(curations$statement_id) ||
               anyNA(curations$statement_id)) {
        "a character statement_id with no NA"
    } else if (!is.numeric(curations$incorrect_count) ||
               anyNA(curations$incorrect_count) ||
               any(curations$incorrect_count < 0) ||
               any(curations$incorrect_count !=
                   round(curations$incorrect_count))) {
        "an incorrect_count of non-negative whole numbers"
    }
    if (!is.null(problem)) {
        stop("get_curations() for ", class(backend), " must return ",
             problem, ".", call. = FALSE)
    }
    invisible(NULL)
}
