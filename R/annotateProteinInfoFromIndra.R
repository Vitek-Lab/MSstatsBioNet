```r
#' Annotate Protein Information from Indra
#'
#' @description This function enriches a data frame with additional protein information sourced from Indra. It appends details such as Uniprot IDs, HGNC IDs, HGNC names, and flags indicating if the protein is a transcription factor, kinase, or phosphatase.
#'
#' @param df \code{data.frame}. The input data frame should be the output of the \code{\link[MSstats]{groupComparison}} function's comparisonResult table. It must contain a list of proteins with their respective p-values, log fold changes (logFCs), and additional columns for HGNC ID and HGNC name.
#' @param proteinIdType \code{character}. Specifies the type of protein ID used in the input data frame. Acceptable values are "Uniprot", "Uniprot_Mnemonic", or "Hgnc_Name".
#' 
#' @return \code{data.frame}. The function returns a data frame with the following columns:
#' \describe{
#'   \item{Protein}{\code{character}. The original protein identifier from the input data frame.}
#'   \item{UniprotID}{\code{character}. The Uniprot ID associated with the protein.}
#'   \item{HgncID}{\code{character}. The HGNC ID associated with the protein.}
#'   \item{HgncName}{\code{character}. The HGNC name associated with the protein.}
#'   \item{IsTranscriptionFactor}{\code{logical}. Indicates whether the protein functions as a transcription factor.}
#'   \item{IsKinase}{\code{logical}. Indicates whether the protein functions as a kinase.}
#'   \item{IsPhosphatase}{\code{logical}. Indicates whether the protein functions as a phosphatase.}
#' }
#' 
#' @examples
#' # Example usage of annotateProteinInfoFromIndra
#' df <- data.frame(Protein = c("CLH1_HUMAN"))
#' annotated_df <- annotateProteinInfoFromIndra(df, "Uniprot_Mnemonic")
#' head(annotated_df)
#' 
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
```