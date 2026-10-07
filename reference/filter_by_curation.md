# Remove evidence curated as incorrect

Subtracts the evidence curated as incorrect, from
[`get_curations()`](https://vitek-lab.github.io/MSstatsBioNet/reference/get_curations.md),
from each edge's `evidence_count`, then drops the edges left below
`min_evidence`, and the nodes those edges leave without any edge. Edges
left with no evidence are dropped whatever `min_evidence` is. A message
says how many edges were dropped.

## Usage

``` r
filter_by_curation(network, min_evidence = 1, backend = NULL)
```

## Arguments

- network:

  list of `nodes` and `edges` from
  [`get_network()`](https://vitek-lab.github.io/MSstatsBioNet/reference/get_network.md).
  `edges` needs the columns `source`, `target`, `statement_id`, and
  `evidence_count`, and `backend_database` unless `backend` is given.

- min_evidence:

  minimum evidence count per edge, after subtracting the incorrect
  evidence

- backend:

  the backend to get the curations from, e.g.
  [`indra_backend()`](https://vitek-lab.github.io/MSstatsBioNet/reference/indra_backend.md).
  `NULL` (default) uses the default backend named in each edge's
  `backend_database`. Pass a backend to use one built with non-default
  settings, such as `indra_backend(curation_url = )`. A backend without
  curations, such as one with no
  [`get_curations()`](https://vitek-lab.github.io/MSstatsBioNet/reference/get_curations.md)
  method, gives an error naming it.

## Value

`network`, with `edges$evidence_count` reduced and the dropped edges and
nodes removed. Other elements of `network` are kept.

## Details

Run it after any filter that pools edges, such as the PTM-site filter of
[`getSubnetworkFromIndra()`](https://vitek-lab.github.io/MSstatsBioNet/reference/getSubnetworkFromIndra.md):
dropping edges first can change which edges such a filter keeps.

## See also

[`get_curations()`](https://vitek-lab.github.io/MSstatsBioNet/reference/get_curations.md)

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
network <- filter_by_curation(network, min_evidence = 2)
#> Dropping 12 edge(s) with fewer than 2 evidence after removing the evidence curated as incorrect.
head(network$edges)
#>    source target interaction directed site confidence evidence_count
#> 4  P05067 O60313     Complex    FALSE <NA>  0.6896123             23
#> 6  P05067 P05362     Complex    FALSE <NA>  0.6549329             13
#> 7  O60313 P05067     Complex    FALSE <NA>  0.6896123             23
#> 9  O00217 O75306     Complex    FALSE <NA>  0.7981153              3
#> 12 O00217 P08574     Complex    FALSE <NA>  0.8216256              3
#> 14 P05362 P05067     Complex    FALSE <NA>  0.6549329             13
#>                                                                evidence_url
#> 4  https://db.indra.bio/statements/from_hash/-12253697417005038?format=html
#> 6  https://db.indra.bio/statements/from_hash/-20220236678417803?format=html
#> 7  https://db.indra.bio/statements/from_hash/-12253697417005038?format=html
#> 9   https://db.indra.bio/statements/from_hash/35453498871697150?format=html
#> 12 https://db.indra.bio/statements/from_hash/-15796466026283356?format=html
#> 14 https://db.indra.bio/statements/from_hash/-20220236678417803?format=html
#>          statement_id backend_database query_type              evidence_sources
#> 4  -12253697417005038            INDRA subnetwork {"sparser": 21, "biogrid": 2}
#> 6  -20220236678417803            INDRA subnetwork   {"sparser": 10, "reach": 3}
#> 7  -12253697417005038            INDRA subnetwork {"sparser": 21, "biogrid": 2}
#> 9   35453498871697150            INDRA subnetwork                {"biogrid": 3}
#> 12 -15796466026283356            INDRA subnetwork                {"biogrid": 3}
#> 14 -20220236678417803            INDRA subnetwork   {"sparser": 10, "reach": 3}
#>    paperCount
#> 4           1
#> 6           1
#> 7           1
#> 9           1
#> 12          1
#> 14          1
# }
```
