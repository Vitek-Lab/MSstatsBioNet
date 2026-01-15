# Populate HGNC IDs in Data Frame

This function populates the HGNC IDs in the data frame based on the
Uniprot IDs.

## Usage

``` r
.populateHgncIdsInDataFrame(df, proteinIdType)
```

## Arguments

- df:

  A data frame containing protein information.

- proteinIdType:

  A character string specifying the type of protein ID. It can be either
  "Uniprot", "Uniprot_Mnemonic", or "Hgnc_Name".

## Value

A data frame with populated HGNC IDs.
