# Copy the values of an entity table onto a network's nodes

Sets `logFC`, `adj.pvalue`, and `node_role` of the nodes that match a
row of `entities`, without querying a backend. Use it to show the
results of another experiment on a network: for example, a network built
from MSstatsPTM's `ADJUSTED.Model`, with its protein nodes colored by
the global profiling results (`PROTEIN.Model`).

## Usage

``` r
add_entity_values(network, entities, overwrite = FALSE)
```

## Arguments

- network:

  list of `nodes` and `edges` from
  [`get_network()`](https://vitek-lab.github.io/MSstatsBioNet/reference/get_network.md)

- entities:

  entity table from
  [`select_entities()`](https://vitek-lab.github.io/MSstatsBioNet/reference/select_entities.md),
  grounded by
  [`convert_ids()`](https://vitek-lab.github.io/MSstatsBioNet/reference/convert_ids.md).
  Its cutoffs set the `node_role` of the matched nodes.

- overwrite:

  whether to replace the `logFC` and `adj.pvalue` of node rows that
  already have them. `FALSE` (default) stops and names those nodes.

## Value

`network`, with the values and roles of `entities` on the matching nodes

## Details

A node row matches an entity row with the same node `id` and `site` (for
a PTM site row, its parent protein and site; for a protein,
`site = NA`). Node rows that match no entity row by `id` are then
matched by grounding (`namespace` and `entity_id`) and `site`. This
finds nodes whose `id` is the backend's name, such as a kinase added
with `get_network(include_entities =)` that INDRA calls `"SRC"`, while
the entity table has `"P12931"`. An entity row that matched a node by
`id` is not used again by grounding, so an isoform node doesn't take the
values of its canonical protein. A node row that matches several entity
rows gets no values, with a message.

A protein that has only PTM site rows in the network, and matches a
protein-level entity row, gets a new protein-level row (`site = NA`)
with that row's values, so
[`cytoscapeNetwork()`](https://vitek-lab.github.io/MSstatsBioNet/reference/cytoscapeNetwork.md)
colors the protein node by it.

Each matched node row gets its `node_role` from its entity row, as set
by
[`select_entities()`](https://vitek-lab.github.io/MSstatsBioNet/reference/select_entities.md):
`"passed_cutoffs"` for a selected row, `"one_condition_missing"` for a
row with an infinite `logFC` selected through `include_infinite_fc`, and
`"user_added"` for a row selected only through `force_include` or not
selected at all. A row that fails the cutoffs is `"user_added"` because
it is shown only because this call put its values on the network. The
new role replaces the node row's old one, so the border and the fill
drawn by
[`cytoscapeNetwork()`](https://vitek-lab.github.io/MSstatsBioNet/reference/cytoscapeNetwork.md)
come from the same row, judged by the cutoffs of the experiment that row
comes from. For example, a protein node filled from the global profiling
results gets a solid border when its global change passes the cutoffs,
while its PTM site nodes keep the roles from the PTM query. Node rows
that match no entity row keep their roles.

Matched node rows get `measured = TRUE`. Rows of `entities` that match
no node are not added:
[`get_network()`](https://vitek-lab.github.io/MSstatsBioNet/reference/get_network.md)
decides which nodes are in the network. Edges and the other elements of
`network` are unchanged. The node table does not record which experiment
the values come from.

## See also

[`filter_by_ptm_site()`](https://vitek-lab.github.io/MSstatsBioNet/reference/filter_by_ptm_site.md),
[`merge_networks()`](https://vitek-lab.github.io/MSstatsBioNet/reference/merge_networks.md)

## Examples

``` r
# A PTM network: phosphorylation of two EGFR sites by SRC, which is not
# in the PTM results, so it has INDRA's name as its id
network <- list(
    nodes = data.frame(
        id = c("P00533", "P00533", "SRC"),
        entity_type = c("ptm_site", "ptm_site", "protein"),
        entity_name = c("EGFR", "EGFR", "SRC"),
        namespace = "HGNC", entity_id = c("3236", "3236", "11283"),
        measured = c(TRUE, TRUE, FALSE), included_in_query = TRUE,
        node_role = c("passed_cutoffs", "passed_cutoffs", "user_added"),
        site = c("S1039", "Y1092", NA),
        has_measured_sites = c(TRUE, TRUE, FALSE),
        logFC = c(1.5, -0.4, NA), adj.pvalue = c(0.01, 0.3, NA)),
    edges = data.frame(
        source = "SRC", target = "P00533", interaction = "Phosphorylation",
        directed = TRUE, site = "Y1092", confidence = 0.9,
        evidence_count = 3L, evidence_url = "https://example.org/1",
        statement_id = "1", backend_database = "INDRA",
        query_type = "subnetwork")
)
# Global profiling results, grounded (convert_ids() would fill these)
global <- prepare_entities(
    data.frame(Protein = c("P00533", "P12931"), log2FC = c(0.2, 1.1),
               adj.pvalue = c(0.6, 0.02)),
    entity_type = "protein", id_type = "uniprot")
global$namespace <- "HGNC"
global$entity_id <- c("3236", "11283")
global$entity_name <- c("EGFR", "SRC")
global <- select_entities(global, pvalue_cutoff = 0.05)
network <- add_entity_values(network, global)
#> Copied the values of entities to 1 node row(s), and added a protein-level row to 1 protein(s) with only PTM site rows.
network$nodes[, c("id", "site", "measured", "node_role", "logFC")]
#>       id  site measured      node_role logFC
#> 1 P00533  <NA>     TRUE     user_added   0.2
#> 2 P00533 S1039     TRUE passed_cutoffs   1.5
#> 3 P00533 Y1092     TRUE passed_cutoffs  -0.4
#> 4    SRC  <NA>     TRUE passed_cutoffs   1.1
```
