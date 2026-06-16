# Populate Transcription Factor Info in Data Frame

Populate Transcription Factor Info in Data Frame

## Usage

``` r
.populateTranscriptionFactorInfoInDataFrame(df, proteinIdType)
```

## Arguments

- df:

  A data frame containing protein information.

- proteinIdType:

  The proteinIdType supplied by the caller. Gene-only flags are `NA` (no
  API call) when this is `"Metabolite"`.

## Value

A data frame with populated transcription factor information.
