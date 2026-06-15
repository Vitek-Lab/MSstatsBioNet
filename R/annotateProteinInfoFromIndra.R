#' Annotate Protein Information from Indra
#'
#' This function standardizes entity identifiers from protein, compound, or
#' gene inputs to a unified namespace using ID conversion from INDRA cogex
#' or Gilda grounding.
#'
#' @param df output of \code{\link[MSstats]{groupComparison}} function's
#'      comparisonResult table. Must contain a \code{Protein} column whose
#'      values are interpreted according to \code{proteinIdType}.
#' @param proteinIdType A character string specifying the type of analyte
#'      identifier in the \code{Protein} column. One of
#'      \code{"Uniprot"}, \code{"Uniprot_Mnemonic"}, \code{"Hgnc_Name"}, or
#'      \code{"Metabolite"}. The \code{"Metabolite"} value treats inputs as
#'      metabolite names and grounds them through Gilda, keeping whatever
#'      namespace Gilda returns (CHEBI / PUBCHEM / CHEMBL / ...).
#' @return A data frame with the following columns:
#' \describe{
#'   \item{Protein}{Character. The original identifier from the input.}
#'   \item{GlobalProtein}{Character. The input identifier with the PTM
#'       site suffix (typically \code{<amino acid><site number>}, e.g.
#'       \code{_S148}) stripped, used as the grounding key.}
#'   \item{UniprotId}{Character. The Uniprot ID of the protein, or
#'       \code{NA} for \code{"Hgnc_Name"} and \code{"Metabolite"} inputs.}
#'   \item{EntityNamespace}{Character. The grounding namespace
#'       (e.g. \code{"HGNC"}, \code{"CHEBI"}). When a single input grounds
#'       to multiple candidates, namespaces are semicolon-joined and
#'       positionally aligned with \code{EntityId} and \code{EntityName}.}
#'   \item{EntityId}{Character. The bare grounding identifier within its
#'       namespace (e.g. \code{"1097"} for HGNC, \code{"28748"} for
#'       CHEBI). Semicolon-joined when multi-grounded.}
#'   \item{EntityName}{Character. The canonical display name from the
#'       grounding source. Semicolon-joined when multi-grounded.}
#'   \item{IsTranscriptionFactor}{Logical. \code{NA} for
#'       \code{proteinIdType == "Metabolite"}.}
#'   \item{IsKinase}{Logical. \code{NA} for
#'       \code{proteinIdType == "Metabolite"}.}
#'   \item{IsPhosphatase}{Logical. \code{NA} for
#'       \code{proteinIdType == "Metabolite"}.}
#' }
#' @examples
#' df <- data.frame(Protein = c("CLH1_HUMAN"))
#' annotated_df <- annotateProteinInfoFromIndra(df, "Uniprot_Mnemonic")
#' head(annotated_df)
#' @export
annotateProteinInfoFromIndra <- function(df, proteinIdType) {
        .validateAnnotateProteinInfoFromIndraInput(df, proteinIdType)
        df <- .populateUniprotIdsInDataFrame(df, proteinIdType)
        df <- .populateEntityInformationInDataFrame(df, proteinIdType)
        df <- .populateTranscriptionFactorInfoInDataFrame(df, proteinIdType)
        df <- .populateKinaseInfoInDataFrame(df, proteinIdType)
        df <- .populatePhophataseInfoInDataFrame(df, proteinIdType)
        return(df)
}

#' Validate Annotate Protein Info Input
#'
#' @param df A data frame containing protein information.
#' @param proteinIdType The proteinIdType supplied by the caller.
#' @return None. Throws an error if validation fails.
.validateAnnotateProteinInfoFromIndraInput <- function(df, proteinIdType) {
        if (!"Protein" %in% colnames(df)) {
                stop("Input dataframe must contain 'Protein' column.")
        }
        allowed <- c("Uniprot", "Uniprot_Mnemonic", "Hgnc_Name", "Metabolite")
        if (length(proteinIdType) != 1 || is.na(proteinIdType) ||
            !(proteinIdType %in% allowed)) {
                stop("Invalid proteinIdType '", proteinIdType, "'. ",
                     "Must be one of: ", paste(allowed, collapse = ", "), ".")
        }
}

#' Populate Uniprot IDs in Data Frame
#'
#' @param df A data frame containing protein information.
#' @param proteinIdType A character string specifying the type of protein ID.
#' @return A data frame with populated Uniprot IDs.
.populateUniprotIdsInDataFrame <- function(df, proteinIdType) {
        if ("GlobalProtein" %in% colnames(df)) {
            protein_ids = unique(as.character(df$GlobalProtein))
        } else {
            df$Protein = as.character(df$Protein)
            df$GlobalProtein = ifelse(grepl("_[A-Z][0-9]", df$Protein),
                                 gsub("_[A-Z][0-9].*", "", df$Protein, perl = TRUE),
                                 df$Protein)
            protein_ids = unique(df$GlobalProtein)
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

        if (proteinIdType == "Hgnc_Name" || proteinIdType == "Metabolite") {
            df$UniprotId <- NA
        }
        return(df)
}

#' Populate Entity Information in Data Frame
#'
#' Initialises the three entity grounding columns and dispatches to the
#' appropriate populator: the INDRA cogex path for UniProt-based inputs,
#' the Gilda grounding path for name-based inputs (HGNC name / metabolite).
#'
#' @param df A data frame containing protein information.
#' @param proteinIdType A character string specifying the type of protein ID.
#' @return A data frame with populated entity grounding columns.
.populateEntityInformationInDataFrame <- function(df, proteinIdType) {
        df$EntityNamespace <- NA
        df$EntityId        <- NA
        df$EntityName      <- NA
        if (proteinIdType == "Uniprot" || proteinIdType == "Uniprot_Mnemonic") {
                df <- .populateEntityInformationWithIndraCogex(df)
        } else {
                df <- .populateEntityInformationWithGilda(df, proteinIdType)
        }
        return(df)
}

#' Populate entity grounding columns via INDRA cogex APIs
#'
#' Converts \code{UniprotId} to an HGNC id via the INDRA cogex endpoint,
#' then looks up the canonical HGNC name. Sets \code{EntityNamespace =
#' "HGNC"} for any row whose UniProt resolved.
#'
#' @param df A data frame with a populated \code{UniprotId} column.
#' @return The data frame with EntityNamespace, EntityId, EntityName set
#'         for resolved rows.
.populateEntityInformationWithIndraCogex <- function(df) {
        validMask <- !is.na(df$UniprotId)
        validUniprots <- unique(df$UniprotId[validMask])
        if (length(validUniprots) > 0) {
                hgncMapping <- .callGetHgncIdsFromUniprotIdsApi(as.list(validUniprots))
                for (uniprotId in names(hgncMapping)) {
                        if (!is.null(hgncMapping[[uniprotId]])) {
                                df$EntityNamespace[df$UniprotId == uniprotId] <- "HGNC"
                                df$EntityId[df$UniprotId == uniprotId]        <- hgncMapping[[uniprotId]]
                        }
                }
        }
        validHgncMask <- !is.na(df$EntityId)
        validHgncs <- unique(df$EntityId[validHgncMask])
        if (length(validHgncs) > 0) {
                nameMapping <- .callGetHgncNamesFromHgncIdsApi(as.list(validHgncs))
                for (entityId in names(nameMapping)) {
                        if (!is.null(nameMapping[[entityId]])) {
                                df$EntityName[df$EntityId == entityId] <- nameMapping[[entityId]]
                        }
                }
        }
        return(df)
}

#' Populate entity grounding columns via Gilda
#'
#' Grounds each \code{GlobalProtein} text through Gilda. For
#' \code{"Hgnc_Name"} the response is filtered to HGNC candidates
#' (and restricted to human via the organism filter); for
#' \code{"Metabolite"} every grounding namespace Gilda returns is kept.
#' Multi-grounded inputs are semicolon-joined and positionally aligned
#' across all three Entity columns.
#'
#' @param df A data frame with a \code{GlobalProtein} column.
#' @param proteinIdType One of \code{"Hgnc_Name"} or \code{"Metabolite"}.
#' @return The data frame with EntityNamespace, EntityId, EntityName set
#'         for resolved rows.
.populateEntityInformationWithGilda <- function(df, proteinIdType) {
        keep_only <- if (proteinIdType == "Hgnc_Name") "HGNC"          else NULL
        organisms <- if (proteinIdType == "Hgnc_Name") list("9606")    else NULL
        textInputs <- unique(df$GlobalProtein)
        if (length(textInputs) > 0) {
                grounding_map <- .callGroundEntitiesFromGildaApi(
                        as.list(textInputs),
                        keep_only = keep_only,
                        organisms = organisms)
                if (!is.null(grounding_map)) {
                        for (input_text in names(grounding_map)) {
                                g <- grounding_map[[input_text]]
                                stopifnot(length(g$ns) == length(g$id),
                                          length(g$ns) == length(g$name))
                                row_mask <- df$GlobalProtein == input_text
                                df$EntityNamespace[row_mask] <- paste(g$ns,   collapse = ";")
                                df$EntityId[row_mask]        <- paste(g$id,   collapse = ";")
                                df$EntityName[row_mask]      <- paste(g$name, collapse = ";")
                        }
                }
        }
        return(df)
}

#' Populate Transcription Factor Info in Data Frame
#'
#' @param df A data frame containing protein information.
#' @param proteinIdType The proteinIdType supplied by the caller. Gene-only
#'        flags are \code{NA} (no API call) when this is \code{"Metabolite"}.
#' @return A data frame with populated transcription factor information.
.populateTranscriptionFactorInfoInDataFrame <- function(df, proteinIdType) {
        df$IsTranscriptionFactor <- NA
        if (proteinIdType == "Metabolite") {
                return(df)
        }
        validNameMask <- !is.na(df$EntityName) & !grepl(";", df$EntityName)
        validNames <- unique(df$EntityName[validNameMask])
        if (length(validNames) > 0) {
                validNamesList <- as.list(validNames)
                charMapping <- .callIsTranscriptionFactorApi(validNamesList)
                for (entityName in names(charMapping)) {
                        if (!is.null(charMapping[[entityName]])) {
                                df$IsTranscriptionFactor[df$EntityName == entityName] <- charMapping[[entityName]]
                        }
                }
        }
        return(df)
}

#' Populate Kinase Info in Data Frame
#'
#' @param df A data frame containing protein information.
#' @param proteinIdType The proteinIdType supplied by the caller. Gene-only
#'        flags are \code{NA} (no API call) when this is \code{"Metabolite"}.
#' @return A data frame with populated kinase information.
.populateKinaseInfoInDataFrame <- function(df, proteinIdType) {
        df$IsKinase <- NA
        if (proteinIdType == "Metabolite") {
                return(df)
        }
        validNameMask <- !is.na(df$EntityName) & !grepl(";", df$EntityName)
        validNames <- unique(df$EntityName[validNameMask])
        if (length(validNames) > 0) {
                validNamesList <- as.list(validNames)
                charMapping <- .callIsKinaseApi(validNamesList)
                for (entityName in names(charMapping)) {
                        if (!is.null(charMapping[[entityName]])) {
                                df$IsKinase[df$EntityName == entityName] <- charMapping[[entityName]]
                        }
                }
        }
        return(df)
}

#' Populate Phosphatase Info in Data Frame
#'
#' @param df A data frame containing protein information.
#' @param proteinIdType The proteinIdType supplied by the caller. Gene-only
#'        flags are \code{NA} (no API call) when this is \code{"Metabolite"}.
#' @return A data frame with populated phosphatase information.
.populatePhophataseInfoInDataFrame <- function(df, proteinIdType) {
        df$IsPhosphatase <- NA
        if (proteinIdType == "Metabolite") {
                return(df)
        }
        validNameMask <- !is.na(df$EntityName) & !grepl(";", df$EntityName)
        validNames <- unique(df$EntityName[validNameMask])
        if (length(validNames) > 0) {
                validNamesList <- as.list(validNames)
                charMapping <- .callIsPhosphataseApi(validNamesList)
                for (entityName in names(charMapping)) {
                        if (!is.null(charMapping[[entityName]])) {
                                df$IsPhosphatase[df$EntityName == entityName] <- charMapping[[entityName]]
                        }
                }
        }
        return(df)
}
