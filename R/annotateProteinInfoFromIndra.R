#' Annotate Protein Information from Indra
#'
#' This function standardizes entity identifiers from protein, compound, or
#' gene inputs to a unified namespace using ID conversion from INDRA cogex
#' or Gilda grounding.
#'
#' @param df output of \code{\link[MSstats]{groupComparison}} function's
#'      comparisonResult table. Must contain a \code{Protein} column whose
#'      values are interpreted according to \code{proteinIdType}. A value
#'      may name a protein group -- several identifiers for the same
#'      quantified analyte joined by \code{";"}, e.g.
#'      \code{"P13747;P23132"} -- in which case every member is grounded
#'      independently and the results are pooled onto the row (see
#'      Details).
#' @param proteinIdType A character string specifying the type of analyte
#'      identifier in the \code{Protein} column. One of
#'      \code{"Uniprot"}, \code{"Uniprot_Mnemonic"}, \code{"Hgnc_Name"}, or
#'      \code{"Metabolite"}. The \code{"Metabolite"} value treats inputs as
#'      metabolite names and grounds them through Gilda, keeping whatever
#'      namespace Gilda returns (CHEBI / PUBCHEM / CHEMBL / ...).
#' @details
#' Protein group members are split on \code{";"}, each member is stripped
#' of its PTM site suffix and grounded on its own, and the groundings of
#' all members are concatenated -- in member order, deduplicated on
#' \code{(EntityNamespace, EntityId)} -- into the semicolon-joined
#' \code{Entity*} columns. This is the same representation used when a
#' single input grounds to several candidates, and
#' \code{\link{getSubnetworkFromIndra}} fans each pair out into its own
#' query node.
#'
#' Because \code{IsTranscriptionFactor} / \code{IsKinase} /
#' \code{IsPhosphatase} describe one gene, they are left \code{NA}
#' whenever a row carries more than one grounding. A group whose members
#' all resolve to the same gene collapses to a single grounding and does
#' get the flags.
#' @return A data frame with the following columns:
#' \describe{
#'   \item{Protein}{Character. The original identifier from the input.}
#'   \item{GlobalProtein}{Character. The input identifier with the PTM
#'       site suffix (typically \code{_<amino acid><site number>}, e.g.
#'       \code{_S148}) stripped from each protein group member, used as
#'       the grounding key. \code{NA} when the input holds no usable
#'       identifier.}
#'   \item{UniprotId}{Character. The Uniprot ID of the protein,
#'       semicolon-joined over the members of a protein group, or
#'       \code{NA} for \code{"Hgnc_Name"} and \code{"Metabolite"} inputs.}
#'   \item{EntityNamespace}{Character. The grounding namespace
#'       (e.g. \code{"HGNC"}, \code{"CHEBI"}). When a row grounds to
#'       multiple candidates -- whether from a protein group or from an
#'       ambiguous single input -- namespaces are semicolon-joined and
#'       positionally aligned with \code{EntityId} and \code{EntityName}.}
#'   \item{EntityId}{Character. The bare grounding identifier within its
#'       namespace (e.g. \code{"1097"} for HGNC, \code{"28748"} for
#'       CHEBI). Semicolon-joined when multi-grounded.}
#'   \item{EntityName}{Character. The canonical display name from the
#'       grounding source. Semicolon-joined when multi-grounded, with
#'       \code{"NA"} in the positions whose name lookup failed.}
#'   \item{IsTranscriptionFactor}{Logical. \code{NA} for
#'       \code{proteinIdType == "Metabolite"} and for multi-grounded rows.}
#'   \item{IsKinase}{Logical. \code{NA} for
#'       \code{proteinIdType == "Metabolite"} and for multi-grounded rows.}
#'   \item{IsPhosphatase}{Logical. \code{NA} for
#'       \code{proteinIdType == "Metabolite"} and for multi-grounded rows.}
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

#' Split a protein group into its member identifiers
#'
#' A \code{Protein} value may name a protein group -- several identifiers
#' for the same quantified analyte joined by \code{";"}. Splits on
#' \code{";"}, trims surrounding whitespace and drops empty members, so a
#' plain single identifier comes back as a length-one vector.
#'
#' @param x A length-one character value, possibly \code{NA}.
#' @return A character vector of member identifiers, empty when the input
#'         holds none.
#' @keywords internal
#' @noRd
.splitProteinGroup <- function(x) {
        if (length(x) == 0 || is.na(x)) {
                return(character(0))
        }
        members <- trimws(unlist(strsplit(as.character(x), ";", fixed = TRUE),
                                 use.names = FALSE))
        return(members[nzchar(members)])
}

#' Join protein group members back into a single value
#'
#' @param members A character vector of member identifiers.
#' @return The members joined by \code{";"}, or \code{NA} when empty.
#' @keywords internal
#' @noRd
.joinProteinGroup <- function(members) {
        if (length(members) == 0) {
                return(NA_character_)
        }
        return(paste(members, collapse = ";"))
}

#' Strip the PTM site suffix from identifiers
#'
#' Removes a trailing \code{_<amino acid><site number>} suffix (e.g.
#' \code{_S148}) from each element, leaving other identifiers untouched.
#'
#' @param x A character vector of identifiers.
#' @return The character vector with site suffixes removed.
#' @keywords internal
#' @noRd
.stripPtmSite <- function(x) {
        return(ifelse(grepl("_[A-Z][0-9]", x),
                      gsub("_[A-Z][0-9].*", "", x, perl = TRUE),
                      x))
}

#' Populate Uniprot IDs in Data Frame
#'
#' Derives \code{GlobalProtein} by stripping the PTM site suffix from each
#' protein group member, then resolves UniProt IDs per member. For a
#' protein group the resolved IDs are semicolon-joined in member order.
#'
#' @param df A data frame containing protein information.
#' @param proteinIdType A character string specifying the type of protein ID.
#' @return A data frame with populated Uniprot IDs.
.populateUniprotIdsInDataFrame <- function(df, proteinIdType) {
        if (!("GlobalProtein" %in% colnames(df))) {
                df$Protein = as.character(df$Protein)
                df$GlobalProtein = vapply(df$Protein, function(protein) {
                        .joinProteinGroup(.stripPtmSite(.splitProteinGroup(protein)))
                }, character(1), USE.NAMES = FALSE)
        }
        df$GlobalProtein = as.character(df$GlobalProtein)
        groupMembers <- lapply(df$GlobalProtein, .splitProteinGroup)
        protein_ids = unique(unlist(groupMembers, use.names = FALSE))
        df$UniprotId <- NA
        if (proteinIdType == "Uniprot") {
                df$UniprotId <- as.character(df$GlobalProtein)
        }

        if (proteinIdType == "Uniprot_Mnemonic") {
                mnemonicProteins <- protein_ids
                if (length(mnemonicProteins) > 0) {
                        uniprotMapping <- .callGetUniprotIdsFromUniprotMnemonicIdsApi(as.list(mnemonicProteins))
                        df$UniprotId <- vapply(groupMembers, function(members) {
                                mapped <- unlist(uniprotMapping[members], use.names = FALSE)
                                .joinProteinGroup(unique(as.character(mapped)))
                        }, character(1), USE.NAMES = FALSE)
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
        df$EntityNamespace <- NA_character_
        df$EntityId        <- NA_character_
        df$EntityName      <- NA_character_
        if (proteinIdType == "Uniprot" || proteinIdType == "Uniprot_Mnemonic") {
                df <- .populateEntityInformationWithIndraCogex(df)
        } else {
                df <- .populateEntityInformationWithGilda(df, proteinIdType)
        }
        return(df)
}

#' Populate entity grounding columns via INDRA cogex APIs
#'
#' Converts each \code{UniprotId} member to an HGNC id via the INDRA cogex
#' endpoint, then looks up the canonical HGNC name. Sets
#' \code{EntityNamespace = "HGNC"} for any row whose UniProt resolved. A
#' protein group contributes one grounding per member that resolved,
#' deduplicated and semicolon-joined across the three Entity columns.
#'
#' @param df A data frame with a populated \code{UniprotId} column, whose
#'        values may be semicolon-joined protein groups.
#' @return The data frame with EntityNamespace, EntityId, EntityName set
#'         for resolved rows.
.populateEntityInformationWithIndraCogex <- function(df) {
        groupMembers <- lapply(df$UniprotId, .splitProteinGroup)
        validUniprots <- unique(unlist(groupMembers, use.names = FALSE))
        if (length(validUniprots) == 0) {
                return(df)
        }
        hgncMapping <- .callGetHgncIdsFromUniprotIdsApi(as.list(validUniprots))
        validHgncs <- unique(as.character(unlist(hgncMapping, use.names = FALSE)))
        nameMapping <- list()
        if (length(validHgncs) > 0) {
                nameResponse <- .callGetHgncNamesFromHgncIdsApi(as.list(validHgncs))
                if (!is.null(nameResponse)) {
                        nameMapping <- nameResponse
                }
        }
        for (i in seq_along(groupMembers)) {
                entityIds <- unique(as.character(
                        unlist(hgncMapping[groupMembers[[i]]], use.names = FALSE)))
                if (length(entityIds) == 0) {
                        next
                }
                entityNames <- vapply(nameMapping[entityIds], function(entityName) {
                        if (is.null(entityName)) NA_character_ else as.character(entityName)[1]
                }, character(1), USE.NAMES = FALSE)
                df$EntityNamespace[i] <- .joinProteinGroup(rep("HGNC", length(entityIds)))
                df$EntityId[i]        <- .joinProteinGroup(entityIds)
                if (!all(is.na(entityNames))) {
                        df$EntityName[i] <- .joinProteinGroup(entityNames)
                }
        }
        return(df)
}

#' Populate entity grounding columns via Gilda
#'
#' Grounds each \code{GlobalProtein} member text through Gilda. For
#' \code{"Hgnc_Name"} the response is filtered to HGNC candidates
#' (and restricted to human via the organism filter); for
#' \code{"Metabolite"} every grounding namespace Gilda returns is kept.
#' Multi-grounded inputs are semicolon-joined and positionally aligned
#' across all three Entity columns; a protein group pools the groundings
#' of all its members into that same representation.
#'
#' @param df A data frame with a \code{GlobalProtein} column, whose values
#'        may be semicolon-joined protein groups.
#' @param proteinIdType One of \code{"Hgnc_Name"} or \code{"Metabolite"}.
#' @return The data frame with EntityNamespace, EntityId, EntityName set
#'         for resolved rows.
.populateEntityInformationWithGilda <- function(df, proteinIdType) {
        keep_only <- if (proteinIdType == "Hgnc_Name") "HGNC"          else NULL
        organisms <- if (proteinIdType == "Hgnc_Name") list("9606")    else NULL
        groupMembers <- lapply(df$GlobalProtein, .splitProteinGroup)
        textInputs <- unique(unlist(groupMembers, use.names = FALSE))
        if (length(textInputs) == 0) {
                return(df)
        }
        grounding_map <- .callGroundEntitiesFromGildaApi(
                as.list(textInputs),
                keep_only = keep_only,
                organisms = organisms)
        if (is.null(grounding_map)) {
                return(df)
        }
        for (i in seq_along(groupMembers)) {
                namespaces <- character(0)
                entityIds <- character(0)
                entityNames <- character(0)
                for (g in grounding_map[groupMembers[[i]]]) {
                        if (is.null(g)) {
                                next
                        }
                        stopifnot(length(g$ns) == length(g$id),
                                  length(g$ns) == length(g$name))
                        namespaces  <- c(namespaces,  as.character(g$ns))
                        entityIds   <- c(entityIds,   as.character(g$id))
                        entityNames <- c(entityNames, as.character(g$name))
                }
                keep <- !duplicated(paste(namespaces, entityIds, sep = ":"))
                if (!any(keep)) {
                        next
                }
                df$EntityNamespace[i] <- .joinProteinGroup(namespaces[keep])
                df$EntityId[i]        <- .joinProteinGroup(entityIds[keep])
                df$EntityName[i]      <- .joinProteinGroup(entityNames[keep])
        }
        return(df)
}

#' Populate Transcription Factor Info in Data Frame
#'
#' Rows carrying more than one grounding -- a semicolon-joined
#' \code{EntityName}, from a protein group or an ambiguous input -- are
#' skipped, because the flag describes a single gene.
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
                                df$IsTranscriptionFactor[which(df$EntityName == entityName)] <- charMapping[[entityName]]
                        }
                }
        }
        return(df)
}

#' Populate Kinase Info in Data Frame
#'
#' Rows carrying more than one grounding -- a semicolon-joined
#' \code{EntityName}, from a protein group or an ambiguous input -- are
#' skipped, because the flag describes a single gene.
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
                                df$IsKinase[which(df$EntityName == entityName)] <- charMapping[[entityName]]
                        }
                }
        }
        return(df)
}

#' Populate Phosphatase Info in Data Frame
#'
#' Rows carrying more than one grounding -- a semicolon-joined
#' \code{EntityName}, from a protein group or an ambiguous input -- are
#' skipped, because the flag describes a single gene.
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
                                df$IsPhosphatase[which(df$EntityName == entityName)] <- charMapping[[entityName]]
                        }
                }
        }
        return(df)
}
