# Ground entities in a backend's namespaces

Fills in `namespace`, `entity_id`, and `entity_name` from each row's
`id` (`parent_id` for `ptm_site` rows), read as the identifier system in
`id_type`. Rows that don't ground are left `NA`.
`backend_capabilities(backend)$id_conversions` lists the (entity type,
identifier system) pairs a backend can convert.

## Usage

``` r
convert_ids(backend, entities, ...)

# S4 method for class 'NetworkBackend'
convert_ids(backend, entities, ...)

# S4 method for class 'IndraBackend'
convert_ids(backend, entities, ...)
```

## Arguments

- backend:

  a `NetworkBackend`, e.g. from
  [`indra_backend()`](https://vitek-lab.github.io/MSstatsBioNet/reference/indra_backend.md)

- entities:

  entity table from
  [`prepare_entities()`](https://vitek-lab.github.io/MSstatsBioNet/reference/prepare_entities.md)

- ...:

  passed to methods

## Value

`entities` with the grounding columns filled in

## Details

Ground the full results table, not only the significant rows, so that
[`get_network()`](https://vitek-lab.github.io/MSstatsBioNet/reference/get_network.md)
can recognize every node that is in the input.

For INDRA, UniProt IDs and mnemonics are mapped through CoGEx, and gene
symbols and chemical names are grounded with Gilda. When an identifier
is a protein group (`"P1;P2"`) or a name with several candidates, the
groundings are `";"`-joined and positionally aligned.

## See also

[`get_entity_properties()`](https://vitek-lab.github.io/MSstatsBioNet/reference/get_entity_properties.md)

## Examples

``` r
df <- data.frame(Protein = c("P04637", "Q00610"))
entities <- prepare_entities(df, entity_type = "protein",
                             id_type = "uniprot")
convert_ids(indra_backend(), entities)
#>       id entity_type id_type namespace entity_id entity_name included_in_query
#> 1 P04637     protein uniprot      HGNC     11998        TP53              TRUE
#> 2 Q00610     protein uniprot      HGNC      2092        CLTC              TRUE
#>   site parent_id organism
#> 1 <NA>      <NA>     9606
#> 2 <NA>      <NA>     9606
```
