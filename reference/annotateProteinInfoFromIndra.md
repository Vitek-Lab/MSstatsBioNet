# Annotate Protein Information from Indra

This function annotates a data frame with protein information from
Indra.

## Usage

``` r
annotateProteinInfoFromIndra(df, proteinIdType)
```

## Arguments

- df:

  output of
  [`groupComparison`](https://rdrr.io/pkg/MSstats/man/groupComparison.html)
  function's comparisonResult table, which contains a list of proteins
  and their corresponding p-values, logFCs, along with additional HGNC
  ID and HGNC name columns

- proteinIdType:

  A character string specifying the type of protein ID. It can be either
  "Uniprot", "Uniprot_Mnemonic", or "Hgnc_Name".

## Value

A data frame with the following columns:

- Protein:

  Character. The original protein identifier.

- UniprotID:

  Character. The Uniprot ID of the protein.

- HgncID:

  Character. The HGNC ID of the protein.

- HgncName:

  Character. The HGNC name of the protein.

- IsTranscriptionFactor:

  Logical. Indicates if the protein is a transcription factor.

- IsKinase:

  Logical. Indicates if the protein is a kinase.

- IsPhosphatase:

  Logical. Indicates if the protein is a phosphatase.

## Examples

``` r
df <- data.frame(Protein = c("CLH1_HUMAN"))
annotated_df <- annotateProteinInfoFromIndra(df, "Uniprot_Mnemonic")
head(annotated_df)
#>      Protein GlobalProtein UniprotId HgncId HgncName IsTranscriptionFactor
#> 1 CLH1_HUMAN    CLH1_HUMAN    Q00610   2092     CLTC                 FALSE
#>   IsKinase IsPhosphatase
#> 1    FALSE         FALSE
```
