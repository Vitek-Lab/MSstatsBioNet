#' Delete an edge from a network edges data frame
#'
#' Removes the row(s) from an edges data frame that match the given
#' \code{source}, \code{target}, and \code{interaction} values.  This is
#' the programmatic counterpart of the interactive Ctrl+click / right-click
#' edge deletion available in \code{\link{cytoscapeNetwork}}.
#'
#' @param edges       Data frame with at minimum columns \code{source},
#'                    \code{target}, and \code{interaction}.
#' @param source      Character. The source node identifier of the edge to
#'                    remove.
#' @param target      Character. The target node identifier of the edge to
#'                    remove.
#' @param interaction Character. The interaction type of the edge to remove.
#'
#' @return The \code{edges} data frame with the matching row(s) removed.
#'
#' @examples
#' edges <- data.frame(
#'   source      = c("TP53",  "MDM2",  "CDKN1A"),
#'   target      = c("MDM2",  "TP53",  "TP53"),
#'   interaction = c("Activation", "Inhibition", "Activation"),
#'   stringsAsFactors = FALSE
#' )
#' deleteEdgeFromNetwork(edges, "MDM2", "TP53", "Inhibition")
#'
#' @export
deleteEdgeFromNetwork <- function(edges, source, target, interaction) {
    if (!is.data.frame(edges)) {
        stop("`edges` must be a data frame.")
    }
    if (!is.character(source) || length(source) != 1L || is.na(source)) {
        stop("`source` must be a single, non-NA character value.")
    }
    if (!is.character(target) || length(target) != 1L || is.na(target)) {
        stop("`target` must be a single, non-NA character value.")
    }
    if (!is.character(interaction) || length(interaction) != 1L || is.na(interaction)) {
        stop("`interaction` must be a single, non-NA character value.")
    }
    required_cols <- c("source", "target", "interaction")
    if (!all(required_cols %in% names(edges))) {
        stop("`edges` must contain columns: source, target, interaction.")
    }
    match_row <- !is.na(edges$source) &
                 !is.na(edges$target) &
                 !is.na(edges$interaction) &
                 edges$source == source &
                 edges$target == target &
                 edges$interaction == interaction
    keep <- !match_row
    edges[keep, , drop = FALSE]
}
}
