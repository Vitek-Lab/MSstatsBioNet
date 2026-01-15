# Populate Uniprot IDs in Data Frame

This function populates the Uniprot IDs in the data frame based on the
protein ID type.

## Usage

``` r
.populateUniprotIdsInDataFrame(df, proteinIdType)
```

## Arguments

- df:

  A data frame containing protein information.

- proteinIdType:

  A character string specifying the type of protein ID. It can be either
  "Uniprot" or "Uniprot_Mnemonic".

## Value

A data frame with populated Uniprot IDs.
