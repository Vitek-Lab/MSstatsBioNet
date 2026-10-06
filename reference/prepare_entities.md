# Prepare an entity table from MSstats results

Builds the table that
[`convert_ids()`](https://vitek-lab.github.io/MSstatsBioNet/reference/convert_ids.md)
grounds,
[`select_entities()`](https://vitek-lab.github.io/MSstatsBioNet/reference/select_entities.md)
flags, and
[`get_network()`](https://vitek-lab.github.io/MSstatsBioNet/reference/get_network.md)
queries. It has one row per analyte, and keeps every analyte, not only
the significant ones, so that nodes returned by a backend can be
recognized as being in the input.

## Usage

``` r
prepare_entities(
  df,
  id_column = "Protein",
  entity_type,
  id_type,
  organism = "9606",
  label = NULL,
  logfc_column = NULL
)
```

## Arguments

- df:

  output of
  [`groupComparison()`](https://rdrr.io/pkg/MSstats/man/groupComparison.html)'s
  `ComparisonResult` table, or any table with one row per analyte.

- id_column:

  name of the column holding the analyte identifiers.

- entity_type:

  one of the entity types (`"protein"`, `"gene"`, `"transcript"`,
  `"ptm_site"`, `"metabolite"`, `"lipid"`, `"drug"`, `"complex"`,
  `"family"`, `"other"`), or the name of a column of `df` holding one
  per row.

- id_type:

  the identifier system of `id_column` (`"uniprot"`,
  `"uniprot_mnemonic"`, `"hgnc_symbol"`, `"ensembl_protein"`,
  `"ensembl_gene"`, `"entrez"`, `"chemical_name"`, `"inchikey"`,
  `"hmdb"`, `"chebi"`, `"chembl"`), or the name of a column of `df`
  holding one per row. Which of them a backend can convert is listed by
  `backend_capabilities(backend)$id_conversions`. For PTM sites, the
  identifier system of the parent protein. `"chemical_name"` is a
  metabolite, lipid, or drug name, common or IUPAC (e.g. `"glucose"`),
  grounded by text matching.

- organism:

  NCBI taxon ID, as a string.

- label:

  the comparison to keep, when `df` has a `Label` column with more than
  one value.

- logfc_column:

  the fold-change column to copy to `logFC`. `NULL` uses whichever of
  `log2FC`, `log10FC`, or `logFC` `df` has. Values are copied unchanged,
  in the log base of the input.

## Value

data.frame with columns `id`, `entity_type`, `id_type`, `namespace`,
`entity_id`, `entity_name` (all `NA` until
[`convert_ids()`](https://vitek-lab.github.io/MSstatsBioNet/reference/convert_ids.md)),
`included_in_query` (`TRUE` until
[`select_entities()`](https://vitek-lab.github.io/MSstatsBioNet/reference/select_entities.md)),
`site` and `parent_id` (for `ptm_site` rows), `organism`, and `logFC`
and `adj.pvalue` when `df` has them. Other columns of `df` are not
copied; join them back on `id`.

## Details

For PTM sites (`entity_type = "ptm_site"`), the site is parsed from the
end of each identifier, e.g. `"P00533_S1039_S1042"` gives parent
`"P00533"` and site `"S1039_S1042"`. A `GlobalProtein` column, as
MSstatsPTM writes, overrides the parsed parent.

## Examples

``` r
input <- data.table::fread(system.file(
    "extdata/groupComparisonModel.csv",
    package = "MSstatsBioNet"
))
entities <- prepare_entities(input, entity_type = "protein",
                             id_type = "uniprot")
head(entities)
#>       id entity_type id_type namespace entity_id entity_name included_in_query
#> 1 O00217     protein uniprot      <NA>      <NA>        <NA>              TRUE
#> 2 O00330     protein uniprot      <NA>      <NA>        <NA>              TRUE
#> 3 O60313     protein uniprot      <NA>      <NA>        <NA>              TRUE
#> 4 O60879     protein uniprot      <NA>      <NA>        <NA>              TRUE
#> 5 O75306     protein uniprot      <NA>      <NA>        <NA>              TRUE
#> 6 P05023     protein uniprot      <NA>      <NA>        <NA>              TRUE
#>   site parent_id organism      logFC  adj.pvalue
#> 1 <NA>      <NA>     9606  2.0285031 0.013821932
#> 2 <NA>      <NA>     9606  1.3000941 0.004457034
#> 3 <NA>      <NA>     9606  0.9299641 0.019584180
#> 4 <NA>      <NA>     9606 -1.9484511 0.001635934
#> 5 <NA>      <NA>     9606  2.4745040 0.004457034
#> 6 <NA>      <NA>     9606  1.8391155 0.003251073

# MSstatsPTM results: the parent protein and site are parsed from the id
ptm <- data.frame(Protein = c("P00533_S1039_S1042", "P04637_S15"),
                  log2FC = c(1.2, -0.8), adj.pvalue = c(0.01, 0.2))
prepare_entities(ptm, entity_type = "ptm_site", id_type = "uniprot")
#>                   id entity_type id_type namespace entity_id entity_name
#> 1 P00533_S1039_S1042    ptm_site uniprot      <NA>      <NA>        <NA>
#> 2         P04637_S15    ptm_site uniprot      <NA>      <NA>        <NA>
#>   included_in_query        site parent_id organism logFC adj.pvalue
#> 1              TRUE S1039_S1042    P00533     9606   1.2       0.01
#> 2              TRUE         S15    P04637     9606  -0.8       0.20
```
