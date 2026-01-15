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
  api_key = ""
)
```

## Arguments

- input:

  output of
  [`groupComparison`](https://rdrr.io/pkg/MSstats/man/groupComparison.html)
  function's comparisionResult table, which contains a list of proteins
  and their corresponding p-values, logFCs, along with additional HGNC
  ID and HGNC name columns

- protein_level_data:

  output of the
  [`dataProcess`](https://rdrr.io/pkg/MSstats/man/dataProcess.html)
  function's ProteinLevelData table, which contains a list of proteins
  and their corresponding abundances. Used for annotating correlation
  information and applying correlation cutoffs.

- pvalueCutoff:

  p-value cutoff for filtering. Default is NULL, i.e. no filtering

- statement_types:

  list of interaction types to filter on. Equivalent to statement type
  in INDRA. Default is NULL.

- paper_count_cutoff:

  number of papers to filter on. Default is 1.

- evidence_count_cutoff:

  number of evidence to filter on for each paper. E.g. A paper may have
  5 sentences describing the same interaction vs 1 sentence. Default is
  1.

- correlation_cutoff:

  if protein_level_abundance is not NULL, apply a cutoff for edges with
  correlation less than a specified cutoff. Default is 0.3

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

- api_key:

  string of INDRA API key for accessing curated statements.

## Value

list of 2 data.frames, nodes and edges

## Examples

``` r
input <- data.table::fread(system.file(
    "extdata/groupComparisonModel.csv",
    package = "MSstatsBioNet"
))
subnetwork <- getSubnetworkFromIndra(input)
#> Warning: NOTICE: This function includes third-party software components
#>         that are licensed under the BSD 2-Clause License. Please ensure to
#>         include the third-party licensing agreements if redistributing this
#>         package or utilizing the results based on this package.
#>         See the LICENSE file for more details.
head(subnetwork$nodes)
#>        id hgncName   Site     logFC  adj.pvalue
#>    <char>   <char> <char>     <num>       <num>
#> 1: O00217   NDUFS8   <NA> 2.0285031 0.013821932
#> 2: O60313     OPA1   <NA> 0.9299641 0.019584180
#> 3: O75306   NDUFS2   <NA> 2.4745040 0.004457034
#> 4: P05023   ATP1A1   <NA> 1.8391155 0.003251073
#> 5: P05067      APP   <NA> 0.7360012 0.020306662
#> 6: P05090     APOD   <NA> 0.5683951 0.013715050
head(subnetwork$edges)
#>   source target site interaction evidenceCount paperCount
#> 1 P05023 O75306 <NA>     Complex             1          1
#> 2 P05067 O60313 <NA>  Activation             2          1
#> 3 O75306 P08574 <NA>     Complex             1          1
#> 4 O60313 O00217 <NA>     Complex             1          1
#> 5 P05362 P05067 <NA>     Complex            13          1
#> 6 O75306 P05067 <NA>     Complex             1          1
#>                                                                                 evidenceLink
#> 1  https://db.indra.bio/statements/from_agents?subject=799@HGNC&object=7708@HGNC&format=html
#> 2  https://db.indra.bio/statements/from_agents?subject=620@HGNC&object=8140@HGNC&format=html
#> 3 https://db.indra.bio/statements/from_agents?subject=7708@HGNC&object=2579@HGNC&format=html
#> 4 https://db.indra.bio/statements/from_agents?subject=8140@HGNC&object=7715@HGNC&format=html
#> 5  https://db.indra.bio/statements/from_agents?subject=5344@HGNC&object=620@HGNC&format=html
#> 6  https://db.indra.bio/statements/from_agents?subject=7708@HGNC&object=620@HGNC&format=html
#>                  sourceCounts          stmt_hash
#> 1              {"biogrid": 1}  -5813063534036006
#> 2                {"reach": 2}   3948742039105656
#> 3              {"biogrid": 1}   6349003830434161
#> 4              {"biogrid": 1} -19747883270157675
#> 5 {"sparser": 10, "reach": 3} -20220236678417803
#> 6              {"biogrid": 1}  22463147519060585
```
