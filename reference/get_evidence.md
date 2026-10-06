# Get the evidence behind network edges

Looks up the evidence sentences and their PubMed IDs for each edge,
using the edge's `statement_id`. Edges that share a `statement_id` get
the same evidence.

## Usage

``` r
get_evidence(backend, edges, ...)

# S4 method for class 'NetworkBackend'
get_evidence(backend, edges, ...)

# S4 method for class 'IndraBackend'
get_evidence(backend, edges, ...)
```

## Arguments

- backend:

  a `NetworkBackend`, e.g. from
  [`indra_backend()`](https://vitek-lab.github.io/MSstatsBioNet/reference/indra_backend.md)

- edges:

  the `edges` of a network from
  [`get_network()`](https://vitek-lab.github.io/MSstatsBioNet/reference/get_network.md).
  Needs the columns `source`, `target`, `interaction`, `site`,
  `evidence_url`, and `statement_id`.

- ...:

  passed to methods

## Value

data.frame with one row per (edge, evidence sentence) pair and the
columns `source`, `target`, `interaction`, `site`, `evidence_url`,
`statement_id`, `text`, and `pmid`. It has no rows, with a warning, when
no edge has evidence text.

## Details

For INDRA, the evidence comes from CoGEx. Evidence without text is left
out, and `pmid` is `""` for evidence from a source with no PubMed ID.

[`filterSubnetworkByContext()`](https://vitek-lab.github.io/MSstatsBioNet/reference/filterSubnetworkByContext.md)
and the topic functions call `get_evidence()` with the backend named in
each edge's `backend_database`, unless their `backend` argument is
given.

## See also

[`filterSubnetworkByContext()`](https://vitek-lab.github.io/MSstatsBioNet/reference/filterSubnetworkByContext.md)

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
evidence <- get_evidence(indra, head(network$edges, 2))
#> Processing 2 unique statement hashes...
#> Fetching evidence for 2 hashes in 1 batch(es) of up to 100...
#> Progress: 1/1 batches (100.0%)
#> Done fetching evidence!
#> Warning: No evidence text found for any statement hash
head(evidence[, c("source", "target", "pmid", "text")])
#> [1] source target pmid   text  
#> <0 rows> (or 0-length row.names)
# }
```
