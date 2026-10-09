# Edge filters shared by every backend's get_network() method.

#' Join the evidence source names of one edge for edges$evidence_sources
#' @param sources character vector of source names
#' @return a single string: the unique names, sorted and \code{";"}-joined,
#' or \code{NA} when there are none
#' @keywords internal
#' @noRd
.join_evidence_sources <- function(sources) {
    sources <- sort(unique(sources[!is.na(sources) & nzchar(sources)]))
    if (length(sources) == 0) NA_character_ else paste(sources, collapse = ";")
}

#' Stop unless min_confidence is NULL or a single number in [0, 1]
#' @param min_confidence the argument of get_network()
#' @keywords internal
#' @noRd
.check_min_confidence <- function(min_confidence) {
    if (is.null(min_confidence)) {
        return(invisible(NULL))
    }
    if (!is.numeric(min_confidence) || length(min_confidence) != 1 ||
        is.na(min_confidence) || min_confidence < 0 || min_confidence > 1) {
        stop("min_confidence must be a single number between 0 and 1.",
             call. = FALSE)
    }
    invisible(NULL)
}

#' Drop edges below a confidence cutoff
#'
#' Edges with no confidence score (\code{NA}) are dropped too, since they
#' can't be shown to pass, and a message says how many.
#' @param edges edges data.frame with a \code{confidence} column
#' @param min_confidence cutoff, or \code{NULL} to keep every edge
#' @return \code{edges}, filtered
#' @keywords internal
#' @noRd
.filter_by_min_confidence <- function(edges, min_confidence) {
    if (is.null(min_confidence)) {
        return(edges)
    }
    no_score <- is.na(edges$confidence)
    if (any(no_score)) {
        message("Dropping ", sum(no_score), " edge(s) with no confidence ",
                "score (NA confidence), because min_confidence is set.")
    }
    edges[!no_score & edges$confidence >= min_confidence, , drop = FALSE]
}

#' Drop edges below an evidence count cutoff
#'
#' Edges with no evidence count (\code{NA}, e.g. from STRING) pass the
#' default \code{min_evidence = 1}, since every edge has at least one piece
#' of evidence. A higher cutoff drops them, since they can't be shown to
#' pass, and a message says how many.
#' @param edges edges data.frame with an \code{evidence_count} column
#' @param min_evidence cutoff
#' @return \code{edges}, filtered
#' @keywords internal
#' @noRd
.filter_by_min_evidence <- function(edges, min_evidence) {
    no_count <- is.na(edges$evidence_count)
    if (min_evidence <= 1) {
        return(edges[no_count | edges$evidence_count >= min_evidence, ,
                     drop = FALSE])
    }
    if (any(no_count)) {
        message("Dropping ", sum(no_count), " edge(s) with no evidence ",
                "count (NA evidence_count), because min_evidence is above 1.")
    }
    edges[!no_count & edges$evidence_count >= min_evidence, , drop = FALSE]
}
