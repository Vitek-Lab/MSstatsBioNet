# Get a network from a backend

Sends the selected entities to a backend and returns the network that
answers the query, as `nodes` and `edges` tables.

## Usage

``` r
get_network(
  backend,
  entities,
  query = subnetwork_query(),
  interaction_types = NULL,
  min_evidence = 1,
  min_confidence = NULL,
  evidence_sources = NULL,
  include_entities = NULL,
  ...
)

# S4 method for class 'NetworkBackend,missing'
get_network(
  backend,
  entities,
  query = subnetwork_query(),
  interaction_types = NULL,
  min_evidence = 1,
  min_confidence = NULL,
  evidence_sources = NULL,
  include_entities = NULL,
  ...
)

# S4 method for class 'NetworkBackend,NetworkQuery'
get_network(
  backend,
  entities,
  query = subnetwork_query(),
  interaction_types = NULL,
  min_evidence = 1,
  min_confidence = NULL,
  evidence_sources = NULL,
  include_entities = NULL,
  ...
)

# S4 method for class 'IndraBackend,SubnetworkQuery'
get_network(
  backend,
  entities,
  query = subnetwork_query(),
  interaction_types = NULL,
  min_evidence = 1,
  min_confidence = NULL,
  evidence_sources = NULL,
  include_entities = NULL,
  ...
)
```

## Arguments

- backend:

  a `NetworkBackend`, e.g. from
  [`indra_backend()`](https://vitek-lab.github.io/MSstatsBioNet/reference/indra_backend.md)

- entities:

  entity table from
  [`prepare_entities()`](https://vitek-lab.github.io/MSstatsBioNet/reference/prepare_entities.md),
  grounded by
  [`convert_ids()`](https://vitek-lab.github.io/MSstatsBioNet/reference/convert_ids.md)
  and flagged by
  [`select_entities()`](https://vitek-lab.github.io/MSstatsBioNet/reference/select_entities.md)

- query:

  a `NetworkQuery` saying which question to ask, e.g.
  [`subnetwork_query()`](https://vitek-lab.github.io/MSstatsBioNet/reference/subnetwork_query.md)
  (the default). See
  [`network_queries`](https://vitek-lab.github.io/MSstatsBioNet/reference/network_queries.md).

- interaction_types:

  values of `edges$interaction` to keep (INDRA statement types, e.g.
  `"Activation"`). `NULL` keeps all.

- min_evidence:

  minimum evidence count per edge

- min_confidence:

  minimum `confidence` per edge, in \[0, 1\]. Edges with no confidence
  score (`NA`) are dropped too, and a message says how many. `NULL`
  applies no cutoff.

- evidence_sources:

  keeps edges with evidence from at least one of these sources, e.g.
  `c("reach")`. `NULL` keeps all.

- include_entities:

  `"namespace:identifier"` groundings to add to the query, e.g.
  `"HGNC:1234"`. Use this for entities outside the input; to keep
  entities of the input that fail the cutoffs, use
  `select_entities(force_include = )`.

- ...:

  passed to methods

## Value

list of `nodes` and `edges` data.frames that meets the contract checked
by
[`validate_network()`](https://vitek-lab.github.io/MSstatsBioNet/reference/validate_network.md)

## Details

Only the rows of `entities` with `included_in_query = TRUE` are sent to
the backend. Every node the backend returns is matched against all rows,
so a node in the input gets `measured = TRUE` and its statistics, and a
node not in the input gets `measured = FALSE` and `NA` statistics.

Each node is one row of `entities` (a PTM site row is drawn on its
parent protein's node). When several rows share a grounding, e.g. two
isoforms that both ground to the same gene, or a protein and a protein
group that contains it, each of their nodes gets the backend's edges,
with its own statistics. So one backend statement can give several
edges, and edges can share a `statement_id`. A statement from a
grounding to itself, such as a homodimer, gives each matching node a
self-loop, and no edges between those nodes. To count statements rather
than edges, count unique `backend_database` and `statement_id` pairs.

`get_network()` prints the question it asks as a message, with the
number of entities, so the query in a saved script or log is readable
without the documentation.

Confidence values are comparable within one backend, not across
backends. For INDRA, `confidence` is the INDRA belief score.

## See also

[`network_queries`](https://vitek-lab.github.io/MSstatsBioNet/reference/network_queries.md),
[`backend_capabilities()`](https://vitek-lab.github.io/MSstatsBioNet/reference/backend_capabilities.md)

## Examples

``` r
input <- data.table::fread(system.file(
    "extdata/groupComparisonModel.csv",
    package = "MSstatsBioNet"
))
entities <- prepare_entities(input, entity_type = "protein",
                             id_type = "uniprot")
# \donttest{
indra <- indra_backend()
entities <- convert_ids(indra, entities)
entities <- select_entities(entities, pvalue_cutoff = 0.05)
network <- get_network(indra, entities, subnetwork_query(),
                       interaction_types = "Complex")
#> INDRA subnetwork: how are 10 selected proteins connected to each other, with no other nodes added?
head(network$nodes)
#>       id entity_type entity_name namespace entity_id measured included_in_query
#> 1 O00217     protein      NDUFS8      HGNC      7715     TRUE              TRUE
#> 2 O60313     protein        OPA1      HGNC      8140     TRUE              TRUE
#> 3 O75306     protein      NDUFS2      HGNC      7708     TRUE              TRUE
#> 4 P05023     protein      ATP1A1      HGNC       799     TRUE              TRUE
#> 5 P05067     protein         APP      HGNC       620     TRUE              TRUE
#> 6 P05090     protein        APOD      HGNC       612     TRUE              TRUE
#>        node_role site has_measured_sites     logFC  adj.pvalue
#> 1 passed_cutoffs <NA>              FALSE 2.0285031 0.013821932
#> 2 passed_cutoffs <NA>              FALSE 0.9299641 0.019584180
#> 3 passed_cutoffs <NA>              FALSE 2.4745040 0.004457034
#> 4 passed_cutoffs <NA>              FALSE 1.8391155 0.003251073
#> 5 passed_cutoffs <NA>              FALSE 0.7360012 0.020306662
#> 6 passed_cutoffs <NA>              FALSE 0.5683951 0.013715050
head(network$edges)
#>   source target interaction directed site confidence evidence_count
#> 1 P05023 O75306     Complex    FALSE <NA>  0.8302067              1
#> 2 O00217 O60313     Complex    FALSE <NA>  0.8302067              1
#> 3 P08574 O75306     Complex    FALSE <NA>  0.8302067              1
#> 4 P05067 O60313     Complex    FALSE <NA>  0.6896123             23
#> 5 O75306 P05067     Complex    FALSE <NA>  0.8302067              1
#> 6 P05067 P05362     Complex    FALSE <NA>  0.6549329             13
#>                                                               evidence_url
#> 1  https://db.indra.bio/statements/from_hash/-5813063534036006?format=html
#> 2 https://db.indra.bio/statements/from_hash/-19747883270157675?format=html
#> 3   https://db.indra.bio/statements/from_hash/6349003830434161?format=html
#> 4 https://db.indra.bio/statements/from_hash/-12253697417005038?format=html
#> 5  https://db.indra.bio/statements/from_hash/22463147519060585?format=html
#> 6 https://db.indra.bio/statements/from_hash/-20220236678417803?format=html
#>         statement_id backend_database query_type              evidence_sources
#> 1  -5813063534036006            INDRA subnetwork                {"biogrid": 1}
#> 2 -19747883270157675            INDRA subnetwork                {"biogrid": 1}
#> 3   6349003830434161            INDRA subnetwork                {"biogrid": 1}
#> 4 -12253697417005038            INDRA subnetwork {"sparser": 21, "biogrid": 2}
#> 5  22463147519060585            INDRA subnetwork                {"biogrid": 1}
#> 6 -20220236678417803            INDRA subnetwork   {"sparser": 10, "reach": 3}
#>   paperCount
#> 1          1
#> 2          1
#> 3          1
#> 4          1
#> 5          1
#> 6          1
# }
```
