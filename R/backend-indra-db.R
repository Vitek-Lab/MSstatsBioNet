# HTTP calls to the INDRA database (db.indra.bio), a different service from
# CoGEx.

#' Count the evidence of a statement curated as incorrect
#'
#' Curations tagged anything other than \code{"correct"} count as
#' incorrect, once per evidence (\code{source_hash}). A failed request
#' warns and counts as 0.
#'
#' @param statement_id INDRA statement hash
#' @param curation_url base URL of the INDRA database
#' @return number of evidences curated as incorrect
#' @keywords internal
#' @noRd
#' @importFrom httr GET status_code content
#' @importFrom jsonlite fromJSON
.get_incorrect_curation_count <- function(statement_id,
                                          curation_url = INDRA_DB_URL) {
    statement_id <- as.character(statement_id)
    url <- file.path(curation_url, "curation/list", statement_id)

    tryCatch({
        response <- GET(url)
        if (status_code(response) == 200) {
            curations <- fromJSON(content(response, "text", encoding = "UTF-8"))
            if (length(curations) == 0) {
                return(0)
            }
            incorrect_curations <- curations[curations$tag != "correct", ]
            length(unique(incorrect_curations$source_hash))
        } else {
            warning(paste("API request failed for hash", statement_id,
                          "with status code", status_code(response)))
            0
        }
    }, error = function(e) {
        warning(paste("Error processing hash", statement_id, ":", e$message))
        0
    })
}
