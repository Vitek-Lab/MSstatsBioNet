# Get the curations of edges from a backend

Curators mark single pieces of evidence for a statement as correct or
incorrect. `get_curations()` counts, for each statement, the evidence
curated as incorrect. INDRA curations come from the INDRA database, with
one request per statement. Most users call
[`filter_by_curation()`](https://vitek-lab.github.io/MSstatsBioNet/reference/filter_by_curation.md),
which subtracts these counts from `evidence_count`.

## Usage

``` r
get_curations(backend, edges, ...)

# S4 method for class 'NetworkBackend'
get_curations(backend, edges, ...)

# S4 method for class 'IndraBackend'
get_curations(backend, edges, ...)
```

## Arguments

- backend:

  a `NetworkBackend`, e.g. from
  [`indra_backend()`](https://vitek-lab.github.io/MSstatsBioNet/reference/indra_backend.md)

- edges:

  the `edges` of a network from
  [`get_network()`](https://vitek-lab.github.io/MSstatsBioNet/reference/get_network.md).
  Needs the column `statement_id`.

- ...:

  passed to methods

## Value

data.frame with one row per unique `statement_id` and the columns
`statement_id` (character) and `incorrect_count` (integer). A failed
request warns and counts as 0. A method for another backend may leave
out statements with no curations, which
[`filter_by_curation()`](https://vitek-lab.github.io/MSstatsBioNet/reference/filter_by_curation.md)
counts as 0, but must return `statement_id` as character:
[`filter_by_curation()`](https://vitek-lab.github.io/MSstatsBioNet/reference/filter_by_curation.md)
errors otherwise, since a numeric hash can lose precision and match no
edge.

## See also

[`filter_by_curation()`](https://vitek-lab.github.io/MSstatsBioNet/reference/filter_by_curation.md)

## Examples

``` r
# \donttest{
input <- data.table::fread(system.file(
    "extdata/groupComparisonModel.csv",
    package = "MSstatsBioNet"
))
indra <- indra_backend()
entities <- prepare_entities(input, entity_type = "protein",
                             id_type = "uniprot")
entities <- convert_ids(indra, entities)
entities <- select_entities(entities, pvalue_cutoff = 0.05)
network <- get_network(indra, entities, interaction_types = "Complex")
#> INDRA subnetwork: how are 10 selected proteins connected to each other, with no other nodes added?
get_curations(indra, head(network$edges, 2))
#>         statement_id incorrect_count
#> 1  -5813063534036006               0
#> 2 -19747883270157675               0
# }
```
