# Keep the edges that modify a PTM site in the network

Keeps an edge when its `site` is one of the PTM sites of its target,
that is, of a node row with the target's `id`. For a network built from
MSstatsPTM results, these are the edges that explain a site in the
input, e.g. a kinase phosphorylating a site whose phosphorylation
changed. Each edge is checked on its own, so the order of this filter
and
[`filter_by_curation()`](https://vitek-lab.github.io/MSstatsBioNet/reference/filter_by_curation.md)
doesn't matter.

## Usage

``` r
filter_by_ptm_site(network)
```

## Arguments

- network:

  list of `nodes` and `edges` from
  [`get_network()`](https://vitek-lab.github.io/MSstatsBioNet/reference/get_network.md),
  with PTM site rows in `nodes` (entity tables with
  `entity_type = "ptm_site"`)

## Value

`network`, with the edges whose site is not a site of their target
removed, and the nodes they leave without edges. Stops when no edges
remain.

## Details

Edges without a `site` are dropped. These include every edge that isn't
a modification (`Complex`, `Activation`, ...), and every edge from a
backend that reports no sites, such as STRING: filter an INDRA network
before merging it with one from such a backend, if those edges should
stay. Nodes are dropped when they had edges and lost them all. Nodes
that had no edge, and other elements of `network`, such as `provenance`,
are kept.

Node sites are split on `"_"` (several sites of one peptide, e.g.
`"S1039_S1042"`) and `";"` (members of a protein group); edge sites on
`";"`.

## See also

[`add_entity_values()`](https://vitek-lab.github.io/MSstatsBioNet/reference/add_entity_values.md)

## Examples

``` r
# A phosphorylation network of CPTAC CCRCC sites and ten kinases, saved
# by inst/script/cptac-ccrcc-network.R
saved <- readRDS(system.file("extdata/cptac-ccrcc-network.rds",
                             package = "MSstatsBioNet"))
network <- filter_by_ptm_site(saved$network)
#> Dropping 580 edge(s) whose site is not a PTM site of their target in the network.
head(network$edges[, c("source", "target", "site")])
#>     source                target site
#> 31   PRKCA  sp|P29353|SHC1_HUMAN S139
#> 144  MAPK1   sp|P29590|PML_HUMAN  S36
#> 147  MAPK1   sp|P29590|PML_HUMAN  S38
#> 148  MAPK1 sp|P78536|ADA17_HUMAN T735
#> 157  MAPK1   sp|P29590|PML_HUMAN  S40
#> 213  PRKCA sp|P78536|ADA17_HUMAN T735
```
