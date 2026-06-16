# Populate Entity Information in Data Frame

Initialises the three entity grounding columns and dispatches to the
appropriate populator: the INDRA cogex path for UniProt-based inputs,
the Gilda grounding path for name-based inputs (HGNC name / metabolite).

## Usage

``` r
.populateEntityInformationInDataFrame(df, proteinIdType)
```

## Arguments

- df:

  A data frame containing protein information.

- proteinIdType:

  A character string specifying the type of protein ID.

## Value

A data frame with populated entity grounding columns.
