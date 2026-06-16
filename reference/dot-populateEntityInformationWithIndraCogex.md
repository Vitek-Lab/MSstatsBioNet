# Populate entity grounding columns via INDRA cogex APIs

Converts `UniprotId` to an HGNC id via the INDRA cogex endpoint, then
looks up the canonical HGNC name. Sets `EntityNamespace = "HGNC"` for
any row whose UniProt resolved.

## Usage

``` r
.populateEntityInformationWithIndraCogex(df)
```

## Arguments

- df:

  A data frame with a populated `UniprotId` column.

## Value

The data frame with EntityNamespace, EntityId, EntityName set for
resolved rows.
