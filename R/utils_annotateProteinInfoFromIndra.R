INDRA_API_URL = "https://discovery.indra.bio"

#' Call API to get UniProt IDs from UniProt mnemonic IDs
#' @param uniprotMnemonicIds list of UniProt mnemonic ids
#' @return list of UniProt IDs
#' @importFrom jsonlite toJSON
#' @importFrom httr POST add_headers content
#' @keywords internal
#' @noRd
.callGetUniprotIdsFromUniprotMnemonicIdsApi <- function(uniprotMnemonicIds) {
    
    if (!is.list(uniprotMnemonicIds)) {
        stop("Input must be a list.")
    }

    if (length(uniprotMnemonicIds) == 0) {
        stop("Input list must not be empty.")
    }
    
    tryCatch({
        # Attempt to convert all elements to character if not already character
        uniprotMnemonicIds <- lapply(uniprotMnemonicIds, function(x) {
            if (!is.character(x)) {
                as.character(x)
            } else {
                x
            }
        })
        
        # Check if conversion was successful
        if (any(!sapply(uniprotMnemonicIds, is.character))) {
            stop("All elements in the list must be character strings representing UniProt mnemonic IDs.")
        }
    }, error = function(e) {
        stop("An error occurred converting uniprot mnemonic IDs to character strings: ", e$message)
    })

    apiUrl <- file.path(INDRA_API_URL, "api/get_uniprot_ids_from_uniprot_mnemonic_ids")

    requestBody <- list(uniprot_mnemonic_ids = uniprotMnemonicIds)
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
    return(res)
}

#' Call API to get HGNC IDs from UniProt IDs
#' @param uniprotIds list of UniProt IDs
#' @return list of HGNC IDs
#' @importFrom jsonlite toJSON
#' @importFrom httr POST add_headers content
#' @keywords internal
#' @noRd
.callGetHgncIdsFromUniprotIdsApi <- function(uniprotIds) {

    if (!is.list(uniprotIds)) {
        stop("Input must be a list.")
    }

    if (any(!sapply(uniprotIds, is.character))) {
        stop("All elements in the list must be character strings representing UniProt IDs.")
    }

    if (length(uniprotIds) == 0) {
        stop("Input list must not be empty.")
    }

    apiUrl <- file.path(INDRA_API_URL, "api/get_hgnc_ids_from_uniprot_ids")

    requestBody <- list(uniprot_ids = uniprotIds)
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
    return(res)
}

#' Call API to get HGNC names from HGNC IDs
#' @param hgncIds list of HGNC IDs
#' @return list of HGNC names
#' @importFrom jsonlite toJSON
#' @importFrom httr POST add_headers content
#' @keywords internal
#' @noRd
.callGetHgncNamesFromHgncIdsApi <- function(hgncIds) {

    if (!is.list(hgncIds)) {
        stop("Input must be a list.")
    }

    if (any(!sapply(hgncIds, is.character))) {
        stop("All elements in the list must be character strings representing HGNC IDs.")
    }

    if (length(hgncIds) == 0) {
        stop("Input list must not be empty.")
    }

    apiUrl <- file.path(INDRA_API_URL, "api/get_hgnc_names_from_hgnc_ids")

    requestBody <- list(hgnc_ids = hgncIds)
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
    return(res)
}

#' Call API to check if genes are kinases
#' @param genes list of gene names
#' @return list indicating if genes are kinases
#' @importFrom jsonlite toJSON
#' @importFrom httr POST add_headers content
#' @keywords internal
#' @noRd
.callIsKinaseApi <- function(genes) {

    if (!is.list(genes)) {
        stop("Input must be a list.")
    }

    if (any(!sapply(genes, is.character))) {
        stop("All elements in the list must be character strings representing gene names.")
    }

    if (length(genes) == 0) {
        stop("Input list must not be empty.")
    }

    apiUrl <- file.path(INDRA_API_URL, "api/is_kinase")

    requestBody <- list(genes = genes)
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
    return(res)
}

#' Call API to check if genes are phosphatases
#' @param genes list of gene names
#' @return list indicating if genes are phosphatases
#' @importFrom jsonlite toJSON
#' @importFrom httr POST add_headers content
#' @keywords internal
#' @noRd
.callIsPhosphataseApi <- function(genes) {

    if (!is.list(genes)) {
        stop("Input must be a list.")
    }

    if (any(!sapply(genes, is.character))) {
        stop("All elements in the list must be character strings representing gene names.")
    }

    if (length(genes) == 0) {
        stop("Input list must not be empty.")
    }

    apiUrl <- file.path(INDRA_API_URL, "api/is_phosphatase")

    requestBody <- list(genes = genes)
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
    return(res)
}

#' Call API to check if genes are transcription factors
#' @param genes list of gene names
#' @return list indicating if genes are transcription factors
#' @importFrom jsonlite toJSON
#' @importFrom httr POST add_headers content
#' @keywords internal
#' @noRd
.callIsTranscriptionFactorApi <- function(genes) {

    if (!is.list(genes)) {
        stop("Input must be a list.")
    }

    if (any(!sapply(genes, is.character))) {
        stop("All elements in the list must be character strings representing gene names.")
    }

    if (length(genes) == 0) {
        stop("Input list must not be empty.")
    }

    apiUrl <- file.path(INDRA_API_URL, "api/is_transcription_factor")

    requestBody <- list(genes = genes)
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
    return(res)
}

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
#' @return Named list keyed by input text. Each value is a list with
#'         three equal-length character vectors: `ns`, `id`, `name`,
#'         positionally aligned across Gilda's returned candidates.
#'         Texts with no surviving grounding are omitted from the result.
#' @importFrom jsonlite toJSON
#' @importFrom httr POST add_headers content
#' @keywords internal
#' @noRd
.callGroundEntitiesFromGildaApi <- function(textInputs, keep_only = NULL) {

    if (!is.list(textInputs)) {
        stop("Input must be a list.")
    }

    if (any(!sapply(textInputs, is.character))) {
        stop("All elements in the list must be character strings.")
    }

    if (length(textInputs) == 0) {
        stop("Input list must not be empty.")
    }

    apiUrl <- file.path("https://grounding.indra.bio/", "ground_multi")

    requestBody <- lapply(textInputs, function(text_input) {
        list(
            text = text_input,
            organisms = list("9606")
        )
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
