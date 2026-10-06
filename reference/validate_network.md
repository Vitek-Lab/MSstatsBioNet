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
have `directed = FALSE`.

Required node column: `id`. An `id` can repeat, once per PTM site row of
the same protein. When present, `entity_type`, `entity_name`,
`namespace`, `entity_id`, `site`, `logFC`, `adj.pvalue`, `measured`,
`included_in_query`, `node_role`, and `has_measured_sites` are
type-checked. Nodes with `measured == FALSE` must have `NA` statistics.

Confidence values are comparable within one `backend_database`, not
across sources.

## Examples

``` r
network <- list(
    nodes = data.frame(id = c("CHK1_HUMAN", "CDC25A_HUMAN")),
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
