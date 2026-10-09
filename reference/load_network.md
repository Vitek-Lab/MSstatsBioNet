# Load a network saved with save_network()

Reads the network and checks it with
[`validate_network()`](https://vitek-lab.github.io/MSstatsBioNet/reference/validate_network.md),
so a file from an older version of the contract gives a clear error here
instead of a failure in a later function.

## Usage

``` r
load_network(file)
```

## Arguments

- file:

  path of an `.rds` file written by
  [`save_network()`](https://vitek-lab.github.io/MSstatsBioNet/reference/save_network.md)

## Value

the network: list of `nodes` and `edges`, plus `provenance` and
`regulators` when it had them

## See also

[`save_network()`](https://vitek-lab.github.io/MSstatsBioNet/reference/save_network.md)

## Examples

``` r
network <- list(
    nodes = data.frame(
        id = c("A", "B"), entity_type = "protein",
        entity_name = c("A", "B"), namespace = "HGNC",
        entity_id = c("1", "2"), measured = TRUE,
        included_in_query = TRUE, node_role = "passed_cutoffs"
    ),
    edges = data.frame(
        source = "A", target = "B", interaction = "Activation",
        directed = TRUE, site = NA_character_, confidence = 0.9,
        evidence_count = 3L, evidence_url = "https://example.org/1",
        statement_id = "1", backend_database = "INDRA",
        query_type = "subnetwork"
    )
)
file <- tempfile(fileext = ".rds")
suppressWarnings(save_network(network, file))
load_network(file)$edges
#>   source target interaction directed site confidence evidence_count
#> 1      A      B  Activation     TRUE <NA>        0.9              3
#>            evidence_url statement_id backend_database query_type
#> 1 https://example.org/1            1            INDRA subnetwork
```
