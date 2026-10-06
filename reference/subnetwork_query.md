# How are my selected entities connected to each other?

Asks which edges connect the selected entities directly, with no other
nodes added. Only entities passed to
[`get_network()`](https://vitek-lab.github.io/MSstatsBioNet/reference/get_network.md)
as `include_entities` are added.

## Usage

``` r
subnetwork_query()
```

## Value

a `SubnetworkQuery` object, to pass to
[`get_network()`](https://vitek-lab.github.io/MSstatsBioNet/reference/get_network.md)

## Details

Uses the rows of `entities` with `included_in_query = TRUE`. Each node
the backend returns is matched against all rows, so that it carries its
statistics. Edges get `query_type = "subnetwork"`. For INDRA, the query
goes to CoGEx `indra_subnetwork_relations`, and takes fewer than 400
groundings (`backend_capabilities(indra_backend())$max_nodes`).

## See also

[`network_queries`](https://vitek-lab.github.io/MSstatsBioNet/reference/network_queries.md)
for the other questions

## Examples

``` r
subnetwork_query()
#> An object of class "SubnetworkQuery"
#> <S4 Type Object>
```
