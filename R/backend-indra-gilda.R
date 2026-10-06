# Gilda grounding call, used by convert_ids(). Moved from
# utils_annotateProteinInfoFromIndra.R in Phase 3 of the API refactor.

#' Call Gilda API to ground entity text against any namespace
#'
#' Posts each input text to Gilda's `ground_multi` endpoint and returns
#' every grounding candidate per input (in Gilda's ranking order). When
#' `keep_only` is set, candidates whose `term$db` does not match are
#' filtered out. The canonical entity name is taken from `term$entry_name`
#' when present, falling back to `term$text` (the input string).
#' @param textInputs list of character strings to ground
#' @param keep_only optional character; if non-NULL, only candidates whose
#'        `term$db == keep_only` are retained
#' @param organisms optional list of NCBI taxonomy ids (e.g.
#'        \code{list("9606")} for human) to constrain Gilda's grounding.
#'        When \code{NULL}, the organisms filter is omitted from the
#'        request body and Gilda may return groundings from any
#'        organism / non-organism namespace (e.g. CHEBI for metabolites).
#' @param grounding_url base URL of Gilda
#' @return Named list keyed by input text. Each value is a list with
#'         three equal-length character vectors: `ns`, `id`, `name`,
#'         positionally aligned across Gilda's returned candidates.
#'         Texts with no surviving grounding are omitted from the result.
#' @importFrom jsonlite toJSON
#' @importFrom httr POST add_headers content
#' @keywords internal
#' @noRd
.callGroundEntitiesFromGildaApi <- function(textInputs, keep_only = NULL,
                                            organisms = NULL,
                                            grounding_url = GILDA_API_URL) {

    if (!is.list(textInputs)) {
        stop("Input must be a list.")
    }

    if (any(!sapply(textInputs, is.character))) {
        stop("All elements in the list must be character strings.")
    }

    if (length(textInputs) == 0) {
        stop("Input list must not be empty.")
    }

    apiUrl <- file.path(grounding_url, "ground_multi")

    requestBody <- lapply(textInputs, function(text_input) {
        entry <- list(text = text_input)
        if (!is.null(organisms)) {
            entry$organisms <- organisms
        }
        entry
    })
    requestBody <- jsonlite::toJSON(requestBody, auto_unbox = TRUE)
    res <- tryCatch({
        response <- POST(
            apiUrl,
            body = requestBody,
            add_headers("Content-Type" = "application/json"),
            encode = "raw"
        )
        content(response)
    }, error = function(e) {
        message("Error in API call: ", e)
        NULL
    })

    if (is.null(res)) {
        return(NULL)
    }

    grounding_map <- list()

    for (i in seq_along(res)) {
        item       <- res[[i]]
        input_text <- as.character(textInputs[[i]])

        ns_vec   <- character(0)
        id_vec   <- character(0)
        name_vec <- character(0)

        for (entry in item) {
            term <- entry$term
            if (is.null(term) || is.null(term$db) || is.null(term$id)) next
            if (!is.null(keep_only) && term$db != keep_only) next

            entry_name <- if (!is.null(term$entry_name) && nzchar(term$entry_name)) {
                term$entry_name
            } else {
                term$text
            }

            ns_vec   <- c(ns_vec,   term$db)
            id_vec   <- c(id_vec,   term$id)
            name_vec <- c(name_vec, entry_name)
        }

        if (length(ns_vec) > 0) {
            grounding_map[[input_text]] <- list(
                ns   = ns_vec,
                id   = id_vec,
                name = name_vec
            )
        }
    }

    return(grounding_map)
}
