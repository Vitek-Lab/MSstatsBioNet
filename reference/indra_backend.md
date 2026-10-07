# Create an INDRA backend

INDRA is a knowledge graph of mechanisms (activations, phosphorylations,
complexes, ...) assembled from the literature and curated databases. The
backend queries INDRA CoGEx for networks and grounds gene symbols and
chemical names with Gilda, INDRA's grounding service. Curations, which
mark evidence as correct or incorrect, come from the INDRA database.
`backend_capabilities(indra_backend())` lists what it supports.

## Usage

``` r
indra_backend(
  cogex_url = INDRA_API_URL,
  grounding_url = GILDA_API_URL,
  curation_url = INDRA_DB_URL
)
```

## Arguments

- cogex_url:

  base URL of INDRA CoGEx

- grounding_url:

  base URL of Gilda

- curation_url:

  base URL of the INDRA database, which holds the curations

## Value

an `IndraBackend` object, to pass to
[`convert_ids()`](https://vitek-lab.github.io/MSstatsBioNet/reference/convert_ids.md),
[`get_entity_properties()`](https://vitek-lab.github.io/MSstatsBioNet/reference/get_entity_properties.md),
[`get_network()`](https://vitek-lab.github.io/MSstatsBioNet/reference/get_network.md),
[`get_evidence()`](https://vitek-lab.github.io/MSstatsBioNet/reference/get_evidence.md),
and
[`get_curations()`](https://vitek-lab.github.io/MSstatsBioNet/reference/get_curations.md)

## Details

This function includes third-party software components that are licensed
under the BSD 2-Clause License. Include the third-party licensing
agreements if redistributing this package or results based on it. See
the LICENSE file for details.

## See also

[`NetworkBackend-class`](https://vitek-lab.github.io/MSstatsBioNet/reference/NetworkBackend-class.md)

## Examples

``` r
indra <- indra_backend()
backend_capabilities(indra)$query_types
#> [1] "subnetwork"
```
