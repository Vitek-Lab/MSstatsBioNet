# Validate a network against the edge and node contract

Checks that a network returned by
[`getSubnetworkFromIndra`](https://vitek-lab.github.io/MSstatsBioNet/reference/getSubnetworkFromIndra.md),
or built by hand from another source, has the columns, types, and
vocabularies that the visualization and filtering functions rely on.

## Usage

``` r
validate_network(network)
```

## Arguments

- network:

  list with `nodes` and `edges` data.frames.

## Value

`network`, invisibly. Stops with an error listing every problem found.

## Details

Required edge columns: `source`, `target` (both matching `nodes$id`),
`interaction` (an INDRA statement type such as `"Activation"` or
`"Complex"`), `directed` (logical), `site` (PTM site on the target such
as `"S148"`, `;`-joined when there are several, or `NA`), `confidence`
(in \[0, 1\], or `NA` when the source provides no score),
`evidence_count` (whole number, at least 1, not `NA`), `evidence_url`,
`statement_id` (character), `backend_database`, and `query_type`. Edges
of the symmetric statement types `"Complex"` and `"Association"` must
have `directed = FALSE`. Several edges can share a `statement_id`: an
undirected statement can be listed in both directions, and a statement
reaches every node with its grounding (see
[`get_network()`](https://vitek-lab.github.io/MSstatsBioNet/reference/get_network.md)).

Required node columns: `id`, `entity_type` (e.g. `"protein"`,
`"ptm_site"`, `"metabolite"`, `"family"`), `entity_name`, `namespace`
and `entity_id` (the grounding, `NA` when unknown), `measured` (logical:
the node is in the input data), `included_in_query` (logical: the node
was part of the query), and `node_role` (why the node is in the network,
e.g. `"passed_cutoffs"` or `"user_added"`, or several joined by `";"`
after
[`merge_networks()`](https://vitek-lab.github.io/MSstatsBioNet/reference/merge_networks.md)).
An `id` can repeat, once per PTM site row of the same protein. When
present, `site`, `has_measured_sites`, `logFC`, and `adj.pvalue` are
type-checked. Nodes with `measured == FALSE` must have `NA` statistics.

Confidence values are comparable within one `backend_database`, not
across sources.

## Examples

``` r
network <- list(
    nodes = data.frame(
        id = c("CHK1_HUMAN", "CDC25A_HUMAN"),
        entity_type = "protein",
        entity_name = c("CHEK1", "CDC25A"),
        namespace = "HGNC",
        entity_id = c("1925", "1725"),
        measured = TRUE,
        included_in_query = TRUE,
        node_role = "passed_cutoffs"
    ),
    edges = data.frame(
        source = "CHK1_HUMAN",
        target = "CDC25A_HUMAN",
        interaction = "Phosphorylation",
        directed = TRUE,
        site = "S76",
        confidence = 0.99,
        evidence_count = 12L,
        evidence_url = paste0("https://db.indra.bio/statements/",
                                "from_hash/-1234?format=html"),
        statement_id = "-1234",
        backend_database = "INDRA",
        query_type = "subnetwork"
    )
)
validate_network(network)
```
