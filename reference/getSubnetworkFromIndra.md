# Get subnetwork from INDRA database

Using differential abundance results from MSstats, this function
retrieves a subnetwork of protein interactions from INDRA database.

## Usage

``` r
getSubnetworkFromIndra(
  input,
  protein_level_data = NULL,
  pvalueCutoff = NULL,
  statement_types = NULL,
  paper_count_cutoff = 1,
  evidence_count_cutoff = 1,
  correlation_cutoff = 0.3,
  sources_filter = NULL,
  logfc_cutoff = NULL,
  force_include_other = NULL,
  filter_by_curation = FALSE,
  filter_by_ptm_site = FALSE,
  include_infinite_fc = FALSE,
  direction = c("both", "up", "down")
)
```

## Arguments

- input:

  output of
  [`groupComparison`](https://rdrr.io/pkg/MSstats/man/groupComparison.html)
  function's comparisionResult table, annotated by
  [`annotateProteinInfoFromIndra`](https://vitek-lab.github.io/MSstatsBioNet/reference/annotateProteinInfoFromIndra.md).
  Must contain `Protein`, `EntityNamespace`, and `EntityId` columns (and
  typically also `EntityName`, `log2FC`, `adj.pvalue`). When an analyte
  grounds to multiple candidates the three `Entity*` columns are
  semicolon-joined and positionally aligned.

- protein_level_data:

  Deprecated, and will be removed in a future release. Output of the
  [`dataProcess`](https://rdrr.io/pkg/MSstats/man/dataProcess.html)
  function's ProteinLevelData table, used to annotate edges with
  correlations and apply `correlation_cutoff`. Supplying it gives a
  deprecation warning.

- pvalueCutoff:

  p-value cutoff for filtering. Default is NULL, i.e. no filtering

- statement_types:

  list of interaction types to filter on. Equivalent to statement type
  in INDRA. Default is NULL.

- paper_count_cutoff:

  Deprecated, and will be removed in a future release. It is ignored:
  paper counts are not available from INDRA, so this filter never had an
  effect for 1 and removed every edge for larger values. Supplying it
  gives a deprecation warning.

- evidence_count_cutoff:

  number of evidence to filter on for each paper. E.g. A paper may have
  5 sentences describing the same interaction vs 1 sentence. Default is
  1.

- correlation_cutoff:

  Deprecated, and will be removed in a future release. If
  `protein_level_data` is not NULL, remove edges whose absolute
  correlation is below this cutoff. Default is 0.3. Supplying it gives a
  deprecation warning.

- sources_filter:

  filtering only on specific sources. Default is no filter, i.e. NULL.
  Otherwise, should be a list, e.g. c('reach', 'medscan').

- logfc_cutoff:

  absolute log fold change cutoff for filtering proteins. Only proteins
  with \|logFC\| greater than this value will be retained. Default is
  NULL, i.e. no logFC filtering.

- force_include_other:

  character vector of identifiers to include in the network, regardless
  if those ids are in the input data. Should be formatted as
  "namespace:identifier", e.g. "HGNC:1234" or "CHEBI:4911".

- filter_by_curation:

  logical, whether to filter out statements that have been curated as
  incorrect in INDRA. Default is FALSE.

- filter_by_ptm_site:

  logical, whether to filter edges based on whether the site information
  from INDRA matches with the PTM site in the input. Default is FALSE.
  Only applicable for differential PTM abundance results.

- include_infinite_fc:

  logical, whether to include proteins with infinite log fold change
  (i.e. proteins that are only detected in one condition). Default is
  FALSE.

- direction:

  Character string specifying the direction of regulation to include.
  One of `"both"` (default), `"up"` (upregulated only), or `"down"`
  (downregulated only).

## Value

list of 2 data.frames, `nodes` and `edges`, that meets the contract
checked by
[`validate_network`](https://vitek-lab.github.io/MSstatsBioNet/reference/validate_network.md).

`edges` has one row per INDRA statement: `source`, `target`,
`interaction` (INDRA statement type), `directed` (`FALSE` for symmetric
types such as `Complex`), `site` (PTM site on the target, or `NA`),
`confidence` (INDRA belief score), `evidence_count`, `evidence_url`
(INDRA page for this statement), `statement_id` (INDRA statement hash,
as character), `backend_database` (`"INDRA"`), `query_type`
(`"subnetwork"`), `evidence_sources` (evidence count per source, as
JSON), and the deprecated `paperCount` (and `correlation` when
`protein_level_data` is given).

`nodes` has one row per analyte (one per site for PTM input, all with
the protein's `id`): `id`, `entity_type` (`"protein"`, `"ptm_site"`,
`"metabolite"`, ...), `entity_name`, `namespace`, `entity_id`,
`measured` (`TRUE` for analytes in `input`), `included_in_query`,
`node_role` (`"passed_cutoffs"`, or `"user_added"` for nodes that are
there only through `force_include_other`), `site`, `has_measured_sites`
(`TRUE` for proteins with PTM site rows), `logFC`, and `adj.pvalue`.
Nodes from `force_include_other` that are not in `input` have
`measured = FALSE` and `NA` statistics.

## Examples

``` r
input <- data.table::fread(system.file(
    "extdata/groupComparisonModel.csv",
    package = "MSstatsBioNet"
))
subnetwork <- getSubnetworkFromIndra(input)
#> INDRA subnetwork: how are 10 selected proteins connected to each other, with no other nodes added?
#> Warning: NOTICE: This function includes third-party software components
#>         that are licensed under the BSD 2-Clause License. Please ensure to
#>         include the third-party licensing agreements if redistributing this
#>         package or utilizing the results based on this package.
#>         See the LICENSE file for more details.
head(subnetwork$nodes)
#>       id entity_type entity_name namespace entity_id measured included_in_query
#> 1 O00217     protein      NDUFS8      HGNC      7715     TRUE              TRUE
#> 2 O60313     protein        OPA1      HGNC      8140     TRUE              TRUE
#> 3 O75306     protein      NDUFS2      HGNC      7708     TRUE              TRUE
#> 4 P05023     protein      ATP1A1      HGNC       799     TRUE              TRUE
#> 5 P05067     protein         APP      HGNC       620     TRUE              TRUE
#> 6 P05090     protein        APOD      HGNC       612     TRUE              TRUE
#>        node_role site has_measured_sites     logFC  adj.pvalue
#> 1 passed_cutoffs <NA>              FALSE 2.0285031 0.013821932
#> 2 passed_cutoffs <NA>              FALSE 0.9299641 0.019584180
#> 3 passed_cutoffs <NA>              FALSE 2.4745040 0.004457034
#> 4 passed_cutoffs <NA>              FALSE 1.8391155 0.003251073
#> 5 passed_cutoffs <NA>              FALSE 0.7360012 0.020306662
#> 6 passed_cutoffs <NA>              FALSE 0.5683951 0.013715050
head(subnetwork$edges)
#>   source target    interaction directed site confidence evidence_count
#> 1 O75306 O00217        Complex    FALSE <NA>  0.7981153              3
#> 2 P08574 O75306        Complex    FALSE <NA>  0.8302067              1
#> 3 P05067 O60313 DecreaseAmount     TRUE <NA>  0.4850742              3
#> 4 P05067 P05362 IncreaseAmount     TRUE <NA>  0.5028536              3
#> 5 P05067 P05090 IncreaseAmount     TRUE <NA>  0.4715316              2
#> 6 P08574 P05023        Complex    FALSE <NA>  0.8302067              1
#>                                                               evidence_url
#> 1  https://db.indra.bio/statements/from_hash/35453498871697150?format=html
#> 2   https://db.indra.bio/statements/from_hash/6349003830434161?format=html
#> 3 https://db.indra.bio/statements/from_hash/-14390698362022978?format=html
#> 4  https://db.indra.bio/statements/from_hash/29503060009374668?format=html
#> 5  https://db.indra.bio/statements/from_hash/13542679335950199?format=html
#> 6 https://db.indra.bio/statements/from_hash/-27854856016676857?format=html
#>         statement_id backend_database query_type evidence_sources paperCount
#> 1  35453498871697150            INDRA subnetwork   {"biogrid": 3}          1
#> 2   6349003830434161            INDRA subnetwork   {"biogrid": 1}          1
#> 3 -14390698362022978            INDRA subnetwork     {"reach": 3}          1
#> 4  29503060009374668            INDRA subnetwork     {"reach": 3}          1
#> 5  13542679335950199            INDRA subnetwork     {"reach": 2}          1
#> 6 -27854856016676857            INDRA subnetwork   {"biogrid": 1}          1
```
