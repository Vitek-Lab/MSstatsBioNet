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

#' Call gilda API to get HGNC IDs from HGNC names
#' @param hgncNames list of hgnc names
#' @return list of HGNC IDs
#' @importFrom jsonlite toJSON
#' @importFrom httr POST add_headers content
#' @keywords internal
#' @noRd
.callGetHgncIdsFromGildaApi <- function(hgncNames) {
    
    if (!is.list(hgncNames)) {
        stop("Input must be a list.")
    }
    
    if (any(!sapply(hgncNames, is.character))) {
        stop("All elements in the list must be character strings representing hgnc names.")
    }
    
    if (length(hgncNames) == 0) {
        stop("Input list must not be empty.")
    }
    
    apiUrl <- file.path("https://grounding.indra.bio/", "ground_multi")
    
    requestBody <- lapply(hgncNames, function(hgnc_name) {
        list(
            text = hgnc_name,
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
    
    hgnc_mapping <- character(0)
    
    for (item in res) {
        # Find the term where db == "HGNC"
        hgnc_term <- NULL
        for (entry in item) {
            if (!is.null(entry$term$db) && entry$term$db == "HGNC") {
                hgnc_term <- entry$term
                break
            }
        }
        
        # Only add to mapping if HGNC term was found
        if (!is.null(hgnc_term)) {
            hgnc_mapping[hgnc_term$text] <- hgnc_term$id
        }
    }
    return(hgnc_mapping)
}
