# Evidence sources of INDRA

The names INDRA gives its evidence sources, in two groups:
`INDRA_DATABASE_SOURCES` for curated databases (e.g. SIGNOR, BioGRID,
PhosphoSitePlus) and `INDRA_TEXT_MINED_SOURCES` for the text-mining
systems that read the literature (e.g. REACH, Sparser). Pass one to
[`get_network`](https://vitek-lab.github.io/MSstatsBioNet/reference/get_network.md)`(evidence_sources = )`
to keep only the edges with evidence from a curated database, or only
those with text-mined evidence.

## Usage

``` r
INDRA_DATABASE_SOURCES

INDRA_TEXT_MINED_SOURCES
```

## Format

character vectors

## Source

`db_sources` and `reader_sources` in INDRA's
`indra.util.statement_presentation` module (INDRA 1.23.0), which INDRA
CoGEx uses for its own source lists. `"bel"` is added to the databases,
as CoGEx does.

## Details

These are the names INDRA uses internally, which are the ones in
`edges$evidence_sources`, and some differ from the full source names:
`"psp"` (PhosphoSitePlus), `"pc"` (Pathway Commons), `"pe"`
(Phospho.ELM), `"vhn"` (VirHostNet), and `"bel_lc"` (BEL large corpus).

## See also

[`get_network()`](https://vitek-lab.github.io/MSstatsBioNet/reference/get_network.md),
[`backend_capabilities()`](https://vitek-lab.github.io/MSstatsBioNet/reference/backend_capabilities.md)

## Examples

``` r
INDRA_DATABASE_SOURCES
#>  [1] "acsn"       "bel"        "bel_lc"     "biogrid"    "cbn"       
#>  [6] "conib"      "creeds"     "crog"       "ctd"        "dgi"       
#> [11] "drugbank"   "hprd"       "minerva"    "omnipath"   "pc"        
#> [16] "pe"         "psp"        "signor"     "tas"        "trrust"    
#> [21] "ubibrowser" "vhn"       
INDRA_TEXT_MINED_SOURCES
#>  [1] "eidos"    "geneways" "gnbr"     "isi"      "medscan"  "reach"   
#>  [7] "rlimsp"   "semrep"   "sparser"  "tees"     "trips"   
# \donttest{
input <- data.table::fread(system.file(
    "extdata/groupComparisonModel.csv",
    package = "MSstatsBioNet"
))
indra <- indra_backend()
entities <- prepare_entities(input, entity_type = "protein",
                             id_type = "uniprot")
entities <- convert_ids(indra, entities)
entities <- select_entities(entities, pvalue_cutoff = 0.05)
curated <- get_network(indra, entities,
                       evidence_sources = INDRA_DATABASE_SOURCES)
#> INDRA subnetwork: how are 10 selected proteins connected to each other, with no other nodes added?
head(curated$edges$evidence_sources)
#> [1] "{\"biogrid\": 1}"                  "{\"biogrid\": 1}"                 
#> [3] "{\"biogrid\": 1}"                  "{\"sparser\": 21, \"biogrid\": 2}"
#> [5] "{\"biogrid\": 1}"                  "{\"sparser\": 21, \"biogrid\": 2}"
# }
```
