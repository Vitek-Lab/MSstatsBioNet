# Populate entity grounding columns via Gilda

Grounds each `GlobalProtein` text through Gilda. For `"Hgnc_Name"` the
response is filtered to HGNC candidates (and restricted to human via the
organism filter); for `"Metabolite"` every grounding namespace Gilda
returns is kept. Multi-grounded inputs are semicolon-joined and
positionally aligned across all three Entity columns.

## Usage

``` r
.populateEntityInformationWithGilda(df, proteinIdType)
```

## Arguments

- df:

  A data frame with a `GlobalProtein` column.

- proteinIdType:

  One of `"Hgnc_Name"` or `"Metabolite"`.

## Value

The data frame with EntityNamespace, EntityId, EntityName set for
resolved rows.
