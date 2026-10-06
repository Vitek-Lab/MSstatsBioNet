# Questions you can ask of a network backend

Each query constructor asks one question of a backend. Pass the query to
[`get_network()`](https://vitek-lab.github.io/MSstatsBioNet/reference/get_network.md).
Every query returns the same `nodes` and `edges` tables (see
[`validate_network()`](https://vitek-lab.github.io/MSstatsBioNet/reference/validate_network.md)),
so visualization and filtering work the same whichever question produced
the network. The `query_type` column of `edges` names the query, without
`_query`.

## Details

|  |  |  |  |
|----|----|----|----|
| **Constructor** | **Question** | **Nodes added** | **`query_type`** |
| [`subnetwork_query()`](https://vitek-lab.github.io/MSstatsBioNet/reference/subnetwork_query.md) | How are my selected entities connected to each other, with no other nodes added? | none, other than `include_entities` | `"subnetwork"` |

More queries (shared regulators, regulator enrichment, paths, and
others) are planned. `backend_capabilities(backend)$query_types` lists
the ones a backend supports.

## Glossary

- entity:

  One row of the entity table from
  [`prepare_entities()`](https://vitek-lab.github.io/MSstatsBioNet/reference/prepare_entities.md):
  one analyte of the input, such as a protein, a PTM site, or a
  metabolite.

- selected:

  An entity with `included_in_query = TRUE`, set by
  [`select_entities()`](https://vitek-lab.github.io/MSstatsBioNet/reference/select_entities.md):
  it passed the cutoffs, or was forced in. Only selected entities are
  sent to the backend.

- in the input (`measured`):

  A node that matches a row of the entity table, whether or not it was
  selected. It carries that row's statistics.

- latent:

  A node that matches no row of the entity table, such as an entity
  added through `include_entities`. It has `measured = FALSE` and `NA`
  statistics. It may still have been measured in another experiment, or
  filtered out before the entity table was built.

- grounding:

  A `"namespace:identifier"` pair that names an entity in a backend,
  such as `"HGNC:11998"` (TP53).

## See also

[`get_network()`](https://vitek-lab.github.io/MSstatsBioNet/reference/get_network.md),
[`backend_capabilities()`](https://vitek-lab.github.io/MSstatsBioNet/reference/backend_capabilities.md)
