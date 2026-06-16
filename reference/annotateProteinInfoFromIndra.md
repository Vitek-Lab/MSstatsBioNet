# Annotate Protein Information from Indra

This function standardizes entity identifiers from protein, compound, or
gene inputs to a unified namespace using ID conversion from INDRA cogex
or Gilda grounding.

## Usage

``` r
annotateProteinInfoFromIndra(df, proteinIdType)
```

## Arguments

- df:

  output of
  [`groupComparison`](https://rdrr.io/pkg/MSstats/man/groupComparison.html)
  function's comparisonResult table. Must contain a `Protein` column
  whose values are interpreted according to `proteinIdType`.

- proteinIdType:

  A character string specifying the type of analyte identifier in the
  `Protein` column. One of `"Uniprot"`, `"Uniprot_Mnemonic"`,
  `"Hgnc_Name"`, or `"Metabolite"`. The `"Metabolite"` value treats
  inputs as metabolite names and grounds them through Gilda, keeping
  whatever namespace Gilda returns (CHEBI / PUBCHEM / CHEMBL / ...).

## Value

A data frame with the following columns:

- Protein:

  Character. The original identifier from the input.

- GlobalProtein:

  Character. The input identifier without the PTM site suffix (typically
  `_<amino acid><site number>`, e.g. `_S148`) stripped, used as the
  grounding key.

- UniprotId:

  Character. The Uniprot ID of the protein, or `NA` for `"Hgnc_Name"`
  and `"Metabolite"` inputs.

- EntityNamespace:

  Character. The grounding namespace (e.g. `"HGNC"`, `"CHEBI"`). When a
  single input grounds to multiple candidates, namespaces are
  semicolon-joined and positionally aligned with `EntityId` and
  `EntityName`.

- EntityId:

  Character. The bare grounding identifier within its namespace (e.g.
  `"1097"` for HGNC, `"28748"` for CHEBI). Semicolon-joined when
  multi-grounded.

- EntityName:

  Character. The canonical display name from the grounding source.
  Semicolon-joined when multi-grounded.

- IsTranscriptionFactor:

  Logical. `NA` for `proteinIdType == "Metabolite"`.

- IsKinase:

  Logical. `NA` for `proteinIdType == "Metabolite"`.

- IsPhosphatase:

  Logical. `NA` for `proteinIdType == "Metabolite"`.

## Examples

``` r
df <- data.frame(Protein = c("CLH1_HUMAN"))
annotated_df <- annotateProteinInfoFromIndra(df, "Uniprot_Mnemonic")
head(annotated_df)
#>      Protein GlobalProtein UniprotId EntityNamespace EntityId EntityName
#> 1 CLH1_HUMAN    CLH1_HUMAN    Q00610            HGNC     2092       CLTC
#>   IsTranscriptionFactor IsKinase IsPhosphatase
#> 1                 FALSE    FALSE         FALSE
```
