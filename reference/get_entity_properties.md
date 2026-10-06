# Add a backend's properties of each entity

Adds one column per property, e.g. `is_kinase`. A property is `NA` for
rows whose `entity_type` it doesn't apply to, and for rows the backend
has no answer for. `backend_capabilities(backend)$entity_properties`
lists the properties a backend supports.

## Usage

``` r
get_entity_properties(backend, entities, properties = NULL, ...)

# S4 method for class 'NetworkBackend'
get_entity_properties(backend, entities, properties = NULL, ...)

# S4 method for class 'IndraBackend'
get_entity_properties(backend, entities, properties = NULL, ...)
```

## Arguments

- backend:

  a `NetworkBackend`, e.g. from
  [`indra_backend()`](https://vitek-lab.github.io/MSstatsBioNet/reference/indra_backend.md)

- entities:

  entity table, grounded by
  [`convert_ids()`](https://vitek-lab.github.io/MSstatsBioNet/reference/convert_ids.md)

- properties:

  the properties to add. `NULL` adds every property the backend
  supports.

- ...:

  passed to methods

## Value

`entities` with one column per property

## Details

For INDRA, the properties are `is_transcription_factor`, `is_kinase`,
and `is_phosphatase`, looked up by gene symbol. Only rows with a single
HGNC grounding get values. A PTM site gets the properties of its parent
protein.

## Examples

``` r
df <- data.frame(Protein = c("P04637", "Q00610"))
indra <- indra_backend()
entities <- prepare_entities(df, entity_type = "protein",
                             id_type = "uniprot")
entities <- convert_ids(indra, entities)
get_entity_properties(indra, entities, properties = "is_kinase")
#>       id entity_type id_type namespace entity_id entity_name included_in_query
#> 1 P04637     protein uniprot      HGNC     11998        TP53              TRUE
#> 2 Q00610     protein uniprot      HGNC      2092        CLTC              TRUE
#>   site parent_id organism is_kinase
#> 1 <NA>      <NA>     9606     FALSE
#> 2 <NA>      <NA>     9606     FALSE
```
