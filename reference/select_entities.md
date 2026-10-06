# Flag the entities to query

Sets `included_in_query` from the statistical cutoffs. Only these rows
are sent to the backend by
[`get_network()`](https://vitek-lab.github.io/MSstatsBioNet/reference/get_network.md).
Rows that fail are kept, so that
[`get_network()`](https://vitek-lab.github.io/MSstatsBioNet/reference/get_network.md)
can still recognize them as being in the input when a backend returns
them. Rows with a missing `adj.pvalue` are not selected.

## Usage

``` r
select_entities(
  entities,
  pvalue_cutoff = NULL,
  logfc_cutoff = NULL,
  direction = c("both", "up", "down"),
  include_infinite_fc = FALSE,
  force_include = NULL
)
```

## Arguments

- entities:

  entity table from
  [`prepare_entities()`](https://vitek-lab.github.io/MSstatsBioNet/reference/prepare_entities.md)

- pvalue_cutoff:

  keep rows with `adj.pvalue` below this. `NULL` applies no cutoff.

- logfc_cutoff:

  keep rows with `abs(logFC)` above this, on the log scale of the input.
  `NULL` applies no cutoff.

- direction:

  `"both"`, `"up"` (`logFC > 0`), or `"down"` (`logFC < 0`).

- include_infinite_fc:

  whether rows with infinite `logFC` (detected in one condition only)
  are selected regardless of `pvalue_cutoff` and `logfc_cutoff`.
  `direction` still applies.

- force_include:

  values of `id`, or `"namespace:identifier"` groundings (e.g.
  `"HGNC:1234"`), selected regardless of the cutoffs. Entities outside
  the table are added with `get_network(include_entities = )` instead.

## Value

`entities` with `included_in_query` set, and `user_added` (`TRUE` for
rows selected only through `force_include`).

## Examples

``` r
input <- data.table::fread(system.file(
    "extdata/groupComparisonModel.csv",
    package = "MSstatsBioNet"
))
entities <- prepare_entities(input, entity_type = "protein",
                             id_type = "uniprot")
entities <- select_entities(entities, pvalue_cutoff = 0.01,
                            direction = "up")
table(entities$included_in_query)
#> 
#> FALSE  TRUE 
#>     7     3 
```
