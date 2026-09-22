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
  whose values are interpreted according to `proteinIdType`. A value may
  name a protein group – several identifiers for the same quantified
  analyte joined by `";"`, e.g. `"P13747;P23132"` – in which case every
  member is grounded independently and the results are pooled onto the
  row.

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
  `_<amino acid><site number>`, e.g. `_S148`) stripped from each protein
  group member, used as the grounding key. `NA` when the input holds no
  usable identifier.

- UniprotId:

  Character. The Uniprot ID of the protein, semicolon-joined in the case
  of multiple proteins, or `NA` for `"Hgnc_Name"` and `"Metabolite"`
  inputs.

- EntityNamespace:

  Character. The grounding namespace (e.g. `"HGNC"`, `"CHEBI"`). When a
  row grounds to multiple candidates – whether from a protein group or
  from an ambiguous single input – namespaces are semicolon-joined and
  positionally aligned with `EntityId` and `EntityName`.

- EntityId:

  Character. The bare grounding identifier within its namespace (e.g.
  `"1097"` for HGNC, `"28748"` for CHEBI). Semicolon-joined when
  multi-grounded.

- EntityName:

  Character. The canonical display name from the grounding source.
  Semicolon-joined when multi-grounded, with `"NA"` in the positions
  whose name lookup failed.

- IsTranscriptionFactor:

  Logical. `NA` for `proteinIdType == "Metabolite"` and for
  multi-grounded rows.

- IsKinase:

  Logical. `NA` for `proteinIdType == "Metabolite"` and for
  multi-grounded rows.

- IsPhosphatase:

  Logical. `NA` for `proteinIdType == "Metabolite"` and for
  multi-grounded rows.

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
