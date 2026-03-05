#' Annotate Protein Information from Indra
#'
#' This function annotates a data frame with protein information from Indra.
#'
#' @param df output of \code{\link[MSstats]{groupComparison}} function's
#'      comparisonResult table, which contains a list of proteins and their
#'      corresponding p-values, logFCs, along with additional HGNC ID and HGNC
#'      name columns
#' @param proteinIdType A character string specifying the type of protein ID.
#'      It can be either "Uniprot", "Uniprot_Mnemonic", or "Hgnc_Name".
#' @return A data frame with the following columns:
#' \describe{
#'   \item{Protein}{Character. The original protein identifier.}
#'   \item{UniprotID}{Character. The Uniprot ID of the protein.}
#'   \item{HgncID}{Character. The HGNC ID of the protein.}
#'   \item{HgncName}{Character. The HGNC name of the protein.}
#'   \item{IsTranscriptionFactor}{Logical. Indicates if the protein is a transcription factor.}
#'   \item{IsKinase}{Logical. Indicates if the protein is a kinase.}
#'   \item{IsPhosphatase}{Logical. Indicates if the protein is a phosphatase.}
#' }
#' @examples
#' df <- data.frame(Protein = c("CLH1_HUMAN"))
#' annotated_df <- annotateProteinInfoFromIndra(df, "Uniprot_Mnemonic")
#' head(annotated_df)
#' @export
annotateProteinInfoFromIndra <- function(df, proteinIdType) {
    .validateAnnotateProteinInfoFromIndraInput(df)
    df <- .populateUniprotIdsInDataFrame(df, proteinIdType)
    df <- .populateHgncIdsInDataFrame(df, proteinIdType)
    df <- .populateHgncNamesInDataFrame(df)
    df <- .populateTranscriptionFactorInfoInDataFrame(df)
    df <- .populateKinaseInfoInDataFrame(df)
    df <- .populatePhophataseInfoInDataFrame(df)
    return(df)
}

#' Validate Annotate Protein Info Input
#'
#' This function validates the input data frame for the annotateProteinInfoFromIndra function.
#'
#' @param df A data frame containing protein information.
#' @return None. Throws an error if validation fails.
.validateAnnotateProteinInfoFromIndraInput <- function(df) {
    if (!"Protein" %in% colnames(df)) {
        stop("Input dataframe must contain 'Protein' column.")
    }
}

#' Populate Uniprot IDs in Data Frame
#'
#' This function populates the Uniprot IDs in the data frame based on the protein ID type.
#'
#' @param df A data frame containing protein information.
#' @param proteinIdType A character string specifying the type of protein ID.
#'        It can be either "Uniprot" or "Uniprot_Mnemonic".
#' @return A data frame with populated Uniprot IDs.
.populateUniprotIdsInDataFrame <- function(df, proteinIdType) {
    if ("GlobalProtein" %in% colnames(df)) {
        protein_ids <- unique(as.character(df$GlobalProtein))
    } else {
        df$Protein <- as.character(df$Protein)
        df$GlobalProtein <- ifelse(grepl("_[A-Z][0-9]", df$Protein),
            gsub("_[A-Z][0-9].*", "", df$Protein, perl = TRUE),
            df$Protein
        )
        protein_ids <- unique(df$GlobalProtein)
    }
    df$UniprotId <- NA
    if (proteinIdType == "Uniprot") {
        df$UniprotId <- as.character(df$GlobalProtein)
    }

    if (proteinIdType == "Uniprot_Mnemonic") {
        mnemonicProteins <- protein_ids
        if (length(mnemonicProteins) > 0) {
            uniprotMapping <- .callGetUniprotIdsFromUniprotMnemonicIdsApi(as.list(mnemonicProteins))
            for (mnemonicId in names(uniprotMapping)) {
                if (!is.null(uniprotMapping[[mnemonicId]])) {
                    df$UniprotId[df$GlobalProtein == mnemonicId] <- uniprotMapping[[mnemonicId]]
                }
            }
        }
    }

    if (proteinIdType == "Hgnc_Name") {
        df$UniprotId <- NA
    }
    return(df)
}

#' Populate HGNC IDs in Data Frame
#'
#' This function populates the HGNC IDs in the data frame based on the Uniprot IDs.
#'
#' @param df A data frame containing protein information.
#' @param proteinIdType A character string specifying the type of protein ID.
#'        It can be either "Uniprot", "Uniprot_Mnemonic", or "Hgnc_Name".
#' @return A data frame with populated HGNC IDs.
.populateHgncIdsInDataFrame <- function(df, proteinIdType) {
    df$HgncId <- NA
    if (proteinIdType == "Uniprot" || proteinIdType == "Uniprot_Mnemonic") {
        validMask <- !is.na(df$UniprotId)
        validUniprots <- unique(df$UniprotId[validMask])
        if (length(validUniprots) > 0) {
            hgncMapping <- .callGetHgncIdsFromUniprotIdsApi(as.list(validUniprots))
            for (uniprotId in names(hgncMapping)) {
                if (!is.null(hgncMapping[[uniprotId]])) {
                    df$HgncId[df$UniprotId == uniprotId] <- hgncMapping[[uniprotId]]
                }
            }
        }
    } else {
        hgncNames <- unique(df$GlobalProtein)
        if (length(hgncNames) > 0) {
            hgncMapping <- .callGetHgncIdsFromGildaApi(as.list(hgncNames))
            for (hgncName in names(hgncMapping)) {
                if (!is.null(hgncMapping[[hgncName]])) {
                    df$HgncId[df$GlobalProtein == hgncName] <- hgncMapping[[hgncName]]
                }
            }
        }
    }

    return(df)
}

#' Populate HGNC Names in Data Frame
#'
#' This function populates the HGNC names in the data frame based on the HGNC IDs.
#'
#' @param df A data frame containing protein information.
#' @return A data frame with populated HGNC names.
.populateHgncNamesInDataFrame <- function(df) {
    df$HgncName <- NA
    validHgncMask <- !is.na(df$HgncId)
    validHgncs <- unique(df$HgncId[validHgncMask])
    if (length(validHgncs) > 0) {
        nameMapping <- .callGetHgncNamesFromHgncIdsApi(as.list(validHgncs))
        for (hgncId in names(nameMapping)) {
            if (!is.null(nameMapping[[hgncId]])) {
                df$HgncName[df$HgncId == hgncId] <- nameMapping[[hgncId]]
            }
        }
    }
    return(df)
}

#' Populate Transcription Factor Info in Data Frame
#'
#' This function populates the transcription factor information in the data frame based on the HGNC names.
#'
#' @param df A data frame containing protein information.
#' @return A data frame with populated transcription factor information.
.populateTranscriptionFactorInfoInDataFrame <- function(df) {
    df$IsTranscriptionFactor <- NA
    validNameMask <- !is.na(df$HgncName)
    validNames <- unique(df$HgncName[validNameMask])
    if (length(validNames) > 0) {
        validNamesList <- as.list(validNames)
        charMapping <- .callIsTranscriptionFactorApi(validNamesList)
        for (hgncName in names(charMapping)) {
            if (!is.null(charMapping[[hgncName]])) {
                df$IsTranscriptionFactor[df$HgncName == hgncName] <- charMapping[[hgncName]]
            }
        }
    }
    return(df)
}

#' Populate Kinase Info in Data Frame
#'
#' This function populates the kinase information in the data frame based on the HGNC names.
#'
#' @param df A data frame containing protein information.
#' @return A data frame with populated kinase information.
.populateKinaseInfoInDataFrame <- function(df) {
    df$IsKinase <- NA
    validNameMask <- !is.na(df$HgncName)
    validNames <- unique(df$HgncName[validNameMask])
    if (length(validNames) > 0) {
        validNamesList <- as.list(validNames)
        charMapping <- .callIsKinaseApi(validNamesList)
        for (hgncName in names(charMapping)) {
            if (!is.null(charMapping[[hgncName]])) {
                df$IsKinase[df$HgncName == hgncName] <- charMapping[[hgncName]]
            }
        }
    }
    return(df)
}

#' Populate Phosphatase Info in Data Frame
#'
#' This function populates the phosphatase information in the data frame based on the HGNC names.
#'
#' @param df A data frame containing protein information.
#' @return A data frame with populated phosphatase information.
.populatePhophataseInfoInDataFrame <- function(df) {
    df$IsPhosphatase <- NA
    validNameMask <- !is.na(df$HgncName)
    validNames <- unique(df$HgncName[validNameMask])
    if (length(validNames) > 0) {
        validNamesList <- as.list(validNames)
        charMapping <- .callIsPhosphataseApi(validNamesList)
        for (hgncName in names(charMapping)) {
            if (!is.null(charMapping[[hgncName]])) {
                df$IsPhosphatase[df$HgncName == hgncName] <- charMapping[[hgncName]]
            }
        }
    }
    return(df)
}
