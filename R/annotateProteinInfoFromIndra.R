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
#'      independently and the results are pooled onto the row.
#' @param proteinIdType A character string specifying the type of analyte
#'      identifier in the \code{Protein} column. One of
#'      \code{"Uniprot"}, \code{"Uniprot_Mnemonic"}, \code{"Hgnc_Name"}, or
#'      \code{"Metabolite"}. The \code{"Metabolite"} value treats inputs as
#'      metabolite names and grounds them through Gilda, keeping whatever
#'      namespace Gilda returns (CHEBI / PUBCHEM / CHEMBL / ...).
#'      
#' @return A data frame with the following columns:
#' \describe{
#'   \item{Protein}{Character. The original identifier from the input.}
#'   \item{GlobalProtein}{Character. The input identifier without the PTM
#'       site suffix (typically \code{_<amino acid><site number>}, e.g.
#'       \code{_S148}) stripped from each protein group member, used as
#'       the grounding key. \code{NA} when the input holds no usable
#'       identifier.}
#'   \item{UniprotId}{Character. The Uniprot ID of the protein,
#'       semicolon-joined in the case of multiple proteins, or
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
        df <- .populateEntityInformationWithIndraBackend(df, proteinIdType)
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

#' Strip the PTM site suffix from identifiers
#'
#' Removes 1+ trailing \code{_<amino acid><site number>} suffixes (e.g.
#' \code{_S148}, \code{_S148_T150}) from each element, 
#' leaving other identifiers untouched.
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
#'
#' @param df A data frame containing protein information.
#' @param proteinIdType A character string specifying the type of protein ID.
#' @return A data frame with populated Uniprot IDs.
#' @noRd
.populateUniprotIdsInDataFrame <- function(df, proteinIdType) {
        if (!("GlobalProtein" %in% colnames(df))) {
                df$Protein = as.character(df$Protein)
                df$GlobalProtein = vapply(df$Protein, function(protein) {
                        .stripPtmSite(protein)
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

#' Populate the grounding and annotation columns through the INDRA backend
#'
#' Builds an entity table with one row per distinct grounding key, runs
#' \code{convert_ids()} and \code{get_annotations()} on it, and copies the
#' results back to every row of \code{df} with that key. The key is
#' \code{UniprotId} for UniProt-based inputs and \code{GlobalProtein}
#' otherwise, so PTM rows are grounded by their stripped parent identifier.
#'
#' @param df A data frame with populated \code{GlobalProtein} and
#'        \code{UniprotId} columns.
#' @param proteinIdType A character string specifying the type of protein ID.
#' @return The data frame with EntityNamespace, EntityId, EntityName,
#'         IsTranscriptionFactor, IsKinase, and IsPhosphatase set.
#' @keywords internal
#' @noRd
.populateEntityInformationWithIndraBackend <- function(df, proteinIdType) {
        usesUniprot <- proteinIdType %in% c("Uniprot", "Uniprot_Mnemonic")
        keys <- as.character(if (usesUniprot) df$UniprotId else df$GlobalProtein)
        hasMembers <- lengths(lapply(keys, .splitProteinGroup)) > 0
        entities <- prepare_entities(
                data.frame(Protein = unique(keys[hasMembers]),
                           stringsAsFactors = FALSE),
                entity_type = if (proteinIdType == "Metabolite") "metabolite" else "protein",
                id_type = switch(proteinIdType,
                                 Uniprot = "uniprot",
                                 Uniprot_Mnemonic = "uniprot",
                                 Hgnc_Name = "hgnc_symbol",
                                 Metabolite = "chemical_name"))
        backend <- indra_backend()
        entities <- convert_ids(backend, entities)
        entities <- get_annotations(backend, entities)
        row <- match(keys, entities$id)
        df$EntityNamespace       <- entities$namespace[row]
        df$EntityId              <- entities$entity_id[row]
        df$EntityName            <- entities$entity_name[row]
        df$IsTranscriptionFactor <- entities$is_transcription_factor[row]
        df$IsKinase              <- entities$is_kinase[row]
        df$IsPhosphatase         <- entities$is_phosphatase[row]
        return(df)
}
