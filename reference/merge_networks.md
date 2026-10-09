# Combine networks from several queries or backends

Combines the networks from several
[`get_network()`](https://vitek-lab.github.io/MSstatsBioNet/reference/get_network.md)
calls, for example two queries against one backend, or one query against
two backends, into one network that meets the same contract
([`validate_network()`](https://vitek-lab.github.io/MSstatsBioNet/reference/validate_network.md)).

## Usage

``` r
merge_networks(..., entities = NULL)
```

## Arguments

- ...:

  networks: lists of `nodes` and `edges` from
  [`get_network()`](https://vitek-lab.github.io/MSstatsBioNet/reference/get_network.md)

- entities:

  `NULL` (default), or the entity table from
  [`prepare_entities()`](https://vitek-lab.github.io/MSstatsBioNet/reference/prepare_entities.md)
  to recompute node status from

## Value

list of `nodes` and `edges`, plus `regulators` and `provenance` when any
network has them

## Details

Edges are the same edge when they have the same `source`, `target`,
`backend_database`, and `statement_id`. One row is kept per edge: the
copy with the largest `evidence_count` (the first one on ties), with the
`query_type` values of all copies joined by `";"`, e.g.
`"subnetwork;mediated"`. A message names the statements whose copies
differ in `evidence_count` or `confidence`. The same relation from two
backends stays two rows, one per `backend_database`, since each has its
own `confidence` and `evidence_url`; confidences are never combined
across backends. Columns that only some networks have are filled with
`NA`.

Nodes are the same node when they have the same `id` and `site` (a
protein has one row per PTM site). Their `node_role` values are joined
by `";"`, e.g. `"passed_cutoffs;mediator"`, `included_in_query` is
`TRUE` if it is in any network, and other columns take the first
non-`NA` value. `has_measured_sites` describes the protein, so it is
`TRUE` on every row of a protein that has it in any network: merging a
protein network with a PTM network marks the protein's protein-level row
too. `included_in_query` describes the query, so it can differ between
networks built from one entity table: a node can be a latent regulator
in one query and part of another.

Without `entities`, `measured`, `logFC`, and `adj.pvalue` must agree for
each node, because networks built from one entity table always agree;
otherwise `merge_networks()` stops and names the nodes. To merge
networks built from different entity tables, pass the entity table as
`entities`: `measured`, `logFC`, `adj.pvalue`, and `has_measured_sites`
are then recomputed from it (nodes not in it are `measured = FALSE` with
`NA` statistics).

The `regulators` and `provenance` tables of the networks, when present,
are combined by row, with columns filled with `NA`. Other list elements
are dropped, with a message.

## See also

[`get_network()`](https://vitek-lab.github.io/MSstatsBioNet/reference/get_network.md),
[`validate_network()`](https://vitek-lab.github.io/MSstatsBioNet/reference/validate_network.md)

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
# Two selections from one entity table: the edges among the proteins
# that pass both cutoffs are in both networks
loose <- get_network(indra, select_entities(entities, pvalue_cutoff = 0.05))
#> INDRA subnetwork: how are 10 selected proteins connected to each other, with no other nodes added?
strict <- get_network(indra, select_entities(entities, pvalue_cutoff = 0.01))
#> INDRA subnetwork: how are 5 selected proteins connected to each other, with no other nodes added?
network <- merge_networks(loose, strict)
# The shared edges come back once
c(loose = nrow(loose$edges), strict = nrow(strict$edges),
  merged = nrow(network$edges))
#>  loose strict merged 
#>     31      2     31 
# }
```
