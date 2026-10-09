# What a backend supports

Lists the queries, identifier conversions, entity properties, and limits
of a backend, so a choice can be checked before a query is sent.

## Usage

``` r
backend_capabilities(backend)

# S4 method for class 'NetworkBackend'
backend_capabilities(backend)

# S4 method for class 'IndraBackend'
backend_capabilities(backend)
```

## Arguments

- backend:

  a `NetworkBackend`, e.g. from
  [`indra_backend()`](https://vitek-lab.github.io/MSstatsBioNet/reference/indra_backend.md)

## Value

named list:

- query_types:

  the `query_type` values of the queries the backend answers, e.g.
  `"subnetwork"` for
  [`subnetwork_query()`](https://vitek-lab.github.io/MSstatsBioNet/reference/subnetwork_query.md)

- id_conversions:

  for each `entity_type`, the `id_type` values
  [`convert_ids()`](https://vitek-lab.github.io/MSstatsBioNet/reference/convert_ids.md)
  can convert

- entity_properties:

  the properties
  [`get_entity_properties()`](https://vitek-lab.github.io/MSstatsBioNet/reference/get_entity_properties.md)
  can add

- interaction_types:

  the `edges$interaction` values the backend returns

- evidence_sources:

  the `evidence_sources` values
  [`get_network()`](https://vitek-lab.github.io/MSstatsBioNet/reference/get_network.md)
  filters on, as a list with `database` and `text_mined` elements, e.g.
  [`INDRA_DATABASE_SOURCES`](https://vitek-lab.github.io/MSstatsBioNet/reference/indra_evidence_sources.md)

- max_nodes:

  for each query type, the largest number of groundings one query can
  take

## See also

[`network_queries`](https://vitek-lab.github.io/MSstatsBioNet/reference/network_queries.md)

## Examples

``` r
backend_capabilities(indra_backend())
#> $query_types
#> [1] "subnetwork"
#> 
#> $id_conversions
#> $id_conversions$protein
#> [1] "uniprot"          "uniprot_mnemonic" "hgnc_symbol"     
#> 
#> $id_conversions$ptm_site
#> [1] "uniprot"          "uniprot_mnemonic" "hgnc_symbol"     
#> 
#> $id_conversions$metabolite
#> [1] "chemical_name"
#> 
#> $id_conversions$lipid
#> [1] "chemical_name"
#> 
#> $id_conversions$drug
#> [1] "chemical_name"
#> 
#> 
#> $entity_properties
#> [1] "is_transcription_factor" "is_kinase"              
#> [3] "is_phosphatase"         
#> 
#> $interaction_types
#>  [1] "Activation"            "Inhibition"            "IncreaseAmount"       
#>  [4] "DecreaseAmount"        "Regulation"            "Influence"            
#>  [7] "Gef"                   "Gap"                   "GtpActivation"        
#> [10] "Modification"          "Phosphorylation"       "Dephosphorylation"    
#> [13] "Autophosphorylation"   "Transphosphorylation"  "Ubiquitination"       
#> [16] "Deubiquitination"      "Sumoylation"           "Desumoylation"        
#> [19] "Hydroxylation"         "Dehydroxylation"       "Acetylation"          
#> [22] "Deacetylation"         "Glycosylation"         "Deglycosylation"      
#> [25] "Farnesylation"         "Defarnesylation"       "Geranylgeranylation"  
#> [28] "Degeranylgeranylation" "Palmitoylation"        "Depalmitoylation"     
#> [31] "Myristoylation"        "Demyristoylation"      "Ribosylation"         
#> [34] "Deribosylation"        "Methylation"           "Demethylation"        
#> [37] "Complex"               "Association"           "Conversion"           
#> [40] "Translocation"        
#> 
#> $evidence_sources
#> $evidence_sources$database
#>  [1] "acsn"       "bel"        "bel_lc"     "biogrid"    "cbn"       
#>  [6] "conib"      "creeds"     "crog"       "ctd"        "dgi"       
#> [11] "drugbank"   "hprd"       "minerva"    "omnipath"   "pc"        
#> [16] "pe"         "psp"        "signor"     "tas"        "trrust"    
#> [21] "ubibrowser" "vhn"       
#> 
#> $evidence_sources$text_mined
#>  [1] "eidos"    "geneways" "gnbr"     "isi"      "medscan"  "reach"   
#>  [7] "rlimsp"   "semrep"   "sparser"  "tees"     "trips"   
#> 
#> 
#> $max_nodes
#> subnetwork 
#>        399 
#> 
```
