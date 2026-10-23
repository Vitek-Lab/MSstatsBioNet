# PTM Analysis with Global Profiling

## Installation

Run this code below to install MSstatsBioNet from bioconductor

``` r

if (!require("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

BiocManager::install("MSstatsBioNet")
```

## Purpose of this vignette

A change in a phosphorylation site means little on its own. The site may
change because its protein’s abundance changed, or because a kinase or
phosphatase acted on it. A global profiling experiment, run on the same
samples, separates the two. This vignette shows how to:

1.  build a network of the kinases known to phosphorylate the sites that
    changed,
2.  keep only the edges that act on a measured site, with
    `filter_by_ptm_site`, and
3.  color each protein by its change in the global profiling experiment,
    with `add_entity_values`.

## Dataset

The data come from the clear cell renal cell carcinoma (CCRCC) cohort of
the Clinical Proteomic Tumor Analysis Consortium (CPTAC), described in
this [paper](https://doi.org/10.1016/j.ccell.2023.06.009) (Li et al.,
Cancer Cell 2023). The study is labeled with TMT and compares tumor
(`T`) with normal adjacent tissue (`N`). It has a
phosphopeptide-enriched experiment and a global profiling experiment on
the same samples. The CPTAC pan-cancer data, including the FragPipe
search results used here, can be downloaded from the [Proteomic Data
Commons](https://pdc.cancer.gov/pdc/cptac-pancancer). The results were
analyzed with MSstatsPTM; the script `inst/script/cptac-ccrcc.R` records
the steps.

MSstatsPTM’s `groupComparisonPTM` returns three tables:

- `PTM.Model`: the change of each site,
- `PROTEIN.Model`: the change of each protein in the global profiling
  experiment,
- `ADJUSTED.Model`: the change of each site adjusted for the change of
  its protein.

We use `ADJUSTED.Model` for the sites. A site that changes only because
its protein does is not a signaling event, and the adjustment removes
it. `PROTEIN.Model` gives the protein-level context.

``` r

ptm = data.table::fread(system.file(
    "extdata/cptac-ccrcc-ptm-adjusted.csv",
    package = "MSstatsBioNet"
))
global = data.table::fread(system.file(
    "extdata/cptac-ccrcc-protein.csv",
    package = "MSstatsBioNet"
))
knitr::kable(head(ptm, 3))
```

| Protein | Label | log2FC | adj.pvalue | GlobalProtein | Adjusted |
|:---|:---|---:|---:|:---|:---|
| sp\|A0MZ66\|SHOT1_HUMAN_S249 | T vs N | 0.59650 | 0.6746 | sp\|A0MZ66\|SHOT1_HUMAN | TRUE |
| sp\|A0MZ66\|SHOT1_HUMAN_S473 | T vs N | -0.41070 | 0.5830 | sp\|A0MZ66\|SHOT1_HUMAN | TRUE |
| sp\|A1IGU5\|ARH37_HUMAN_S15 | T vs N | 0.04957 | 0.9747 | sp\|A1IGU5\|ARH37_HUMAN | TRUE |

``` r

knitr::kable(head(global, 3))
```

| Protein                     | Label  |  log2FC | adj.pvalue |
|:----------------------------|:-------|--------:|-----------:|
| sp\|A0A024RBG1\|NUD4B_HUMAN | T vs N | -0.8327 |     0.1436 |
| sp\|A0A075B6K0\|LV316_HUMAN | T vs N |  1.1660 |     0.2180 |
| sp\|A0A075B6K5\|LV39_HUMAN  | T vs N |  1.5630 |     0.5190 |

The identifiers are in FragPipe’s form,
`sp|<UniProt accession>|<UniProt mnemonic>`, followed by the site for
PTMs. `GlobalProtein` names the protein of each site. Sites whose
protein was not measured in the global profiling experiment cannot be
adjusted: MSstatsPTM sets their `GlobalProtein` to `"missing"` and
`Adjusted` to `FALSE`, and reports their unadjusted change.

``` r

table(Adjusted = ptm$Adjusted)
#> Adjusted
#> FALSE  TRUE 
#>  1141  1972
```

## Building the entity tables

`prepare_entities` builds one entity table per experiment. For the
sites, `entity_type = "ptm_site"` parses the site from the end of each
identifier and takes the parent protein from `GlobalProtein` (or from
the identifier, where it is `"missing"`). FragPipe identifiers can be
used as they are with `id_type = "uniprot"`.

``` r

library(MSstatsBioNet)
#> Loading required package: MSstats
#> 
#> Attaching package: 'MSstats'
#> The following object is masked from 'package:grDevices':
#> 
#>     savePlot
sites <- prepare_entities(ptm, entity_type = "ptm_site", id_type = "uniprot")
proteins <- prepare_entities(global, entity_type = "protein",
                             id_type = "uniprot")
knitr::kable(head(sites[, c("id", "parent_id", "site", "logFC",
                            "adj.pvalue")], 3))
```

| id | parent_id | site | logFC | adj.pvalue |
|:---|:---|:---|---:|---:|
| sp\|A0MZ66\|SHOT1_HUMAN_S249 | sp\|A0MZ66\|SHOT1_HUMAN | S249 | 0.59650 | 0.6746 |
| sp\|A0MZ66\|SHOT1_HUMAN_S473 | sp\|A0MZ66\|SHOT1_HUMAN | S473 | -0.41070 | 0.5830 |
| sp\|A1IGU5\|ARH37_HUMAN_S15 | sp\|A1IGU5\|ARH37_HUMAN | S15 | 0.04957 | 0.9747 |

## ID Conversion

`convert_ids` maps each UniProt accession to its HGNC gene, which INDRA
uses. Sites are mapped through their parent protein.

So that this vignette builds without network access, the chunks that
call INDRA (`convert_ids` and `get_network`) are not run when the
vignette is built. Their results were saved with
`inst/script/cptac-ccrcc-network.R` and are loaded from
`extdata/cptac-ccrcc-network.rds` instead. The output below comes from
those saved results.

``` r

indra <- indra_backend()
sites <- convert_ids(indra, sites)
proteins <- convert_ids(indra, proteins)
```

``` r

knitr::kable(head(sites[, c("id", "namespace", "entity_id",
                            "entity_name")], 3))
```

| id                           | namespace | entity_id | entity_name |
|:-----------------------------|:----------|:----------|:------------|
| sp\|A0MZ66\|SHOT1_HUMAN_S249 | HGNC      | 29319     | SHTN1       |
| sp\|A0MZ66\|SHOT1_HUMAN_S473 | HGNC      | 29319     | SHTN1       |
| sp\|A1IGU5\|ARH37_HUMAN_S15  | HGNC      | 34430     | ARHGEF37    |

``` r

c(sites = sum(!is.na(sites$entity_id)), of = nrow(sites))
#> sites    of 
#>  3113  3113
c(proteins = sum(!is.na(proteins$entity_id)), of = nrow(proteins))
#> proteins       of 
#>     4404     4405
```

## Selecting Sites

`select_entities` flags the sites to send to INDRA. Here we keep the
sites with an adjusted p-value below 0.05.

``` r

selected <- select_entities(sites, pvalue_cutoff = 0.05)
table(selected$included_in_query)
#> 
#> FALSE  TRUE 
#>  2877   236
```

## Network Query

### Sites only

`get_network` with
[`subnetwork_query()`](https://vitek-lab.github.io/MSstatsBioNet/reference/subnetwork_query.md)
asks how the proteins of the selected sites are connected to each other.
We keep only `Phosphorylation` statements.

``` r

network_sites <- get_network(indra, selected, subnetwork_query(),
                             interaction_types = "Phosphorylation")
#> INDRA subnetwork: how are 236 selected PTM sites connected to each other, with no other nodes added?
```

``` r

nrow(network_sites$edges)
#> [1] 30
```

INDRA reports the site each phosphorylation acts on, in the `site`
column of the edges. Most of these edges act on a site that was not
measured, or that did not change. `filter_by_ptm_site` keeps only the
edges whose site is one of the sites of their target in the network.

``` r

filtered_sites <- filter_by_ptm_site(network_sites)
#> Dropping 29 edge(s) whose site is not a PTM site of their target in the network.
knitr::kable(filtered_sites$edges[, c("source", "target", "interaction",
                                      "site", "evidence_count")])
```

|  | source | target | interaction | site | evidence_count |
|:---|:---|:---|:---|:---|---:|
| 21 | sp\|Q07002\|CDK18_HUMAN | sp\|P41181\|AQP2_HUMAN | Phosphorylation | S261 | 4 |

Few of the proteins whose sites changed are known to phosphorylate each
other, so one edge remains.

### Adding kinases with `include_entities`

The kinases that act on these sites may not have changed sites
themselves, so they are not selected. They can be added to the query
with `include_entities`, using `"namespace:identifier"` IDs. Here we add
ten kinases that act on many substrates.

``` r

kinases <- c(MTOR = "HGNC:3942", AKT1 = "HGNC:391", MAPK1 = "HGNC:6871",
             MAPK3 = "HGNC:6877", CDK1 = "HGNC:1722", CSNK2A1 = "HGNC:2457",
             PRKACA = "HGNC:9380", GSK3B = "HGNC:4617", SRC = "HGNC:11283",
             PRKCA = "HGNC:9393")
```

``` r

network <- get_network(indra, selected, subnetwork_query(),
                       interaction_types = "Phosphorylation",
                       include_entities = unname(kinases))
#> INDRA subnetwork: how are 236 selected PTM sites and 10 added entities connected to each other, with no other nodes added?
```

``` r

nrow(network$edges)
#> [1] 599
network <- filter_by_ptm_site(network)
#> Dropping 580 edge(s) whose site is not a PTM site of their target in the network.
knitr::kable(network$edges[, c("source", "target", "site",
                               "evidence_count")])
```

|     | source                  | target                  | site  | evidence_count |
|:----|:------------------------|:------------------------|:------|---------------:|
| 31  | PRKCA                   | sp\|P29353\|SHC1_HUMAN  | S139  |             10 |
| 144 | MAPK1                   | sp\|P29590\|PML_HUMAN   | S36   |              3 |
| 147 | MAPK1                   | sp\|P29590\|PML_HUMAN   | S38   |              2 |
| 148 | MAPK1                   | sp\|P78536\|ADA17_HUMAN | T735  |              3 |
| 157 | MAPK1                   | sp\|P29590\|PML_HUMAN   | S40   |              2 |
| 213 | PRKCA                   | sp\|P78536\|ADA17_HUMAN | T735  |              1 |
| 235 | AKT1                    | sp\|P02545\|LMNA_HUMAN  | S301  |              1 |
| 249 | MAPK3                   | sp\|P29590\|PML_HUMAN   | S36   |              2 |
| 254 | MAPK3                   | sp\|P29590\|PML_HUMAN   | S38   |              1 |
| 257 | PRKCA                   | sp\|P41181\|AQP2_HUMAN  | S256  |              1 |
| 266 | MAPK3                   | sp\|P29590\|PML_HUMAN   | S40   |              1 |
| 284 | PRKACA                  | sp\|Q16665\|HIF1A_HUMAN | T455  |              1 |
| 316 | MAPK3                   | sp\|P78536\|ADA17_HUMAN | T735  |              6 |
| 332 | CDK1                    | sp\|P41181\|AQP2_HUMAN  | S261  |              1 |
| 345 | CDK1                    | sp\|O60885\|BRD4_HUMAN  | S1117 |              2 |
| 422 | CSNK2A1                 | sp\|P25788\|PSA3_HUMAN  | S250  |              1 |
| 462 | SRC                     | sp\|P46940\|IQGA1_HUMAN | Y1510 |              2 |
| 505 | sp\|Q07002\|CDK18_HUMAN | sp\|P41181\|AQP2_HUMAN  | S261  |              4 |
| 513 | PRKACA                  | sp\|P41181\|AQP2_HUMAN  | S256  |              2 |

The added kinases are not in the PTM entity table, so their nodes have
INDRA’s gene symbol as `id`, `measured = FALSE`, and no statistics.

Run `filter_by_ptm_site` before `merge_networks` when combining with a
backend that reports no sites: it drops every edge without a site.

## Adding global profiling values

`add_entity_values` copies the `logFC` and `adj.pvalue` of an entity
table onto the matching nodes of a network, without a new query. With
the global profiling table:

- a protein with only site rows in the network gets a protein-level row,
  with the protein’s change;
- a kinase added with `include_entities` is matched by its HGNC
  grounding, since its `id` is INDRA’s gene symbol, and gets its change
  too.

The entity table must come from `convert_ids`, so that nodes can be
matched by grounding.

``` r

network <- add_entity_values(network, proteins)
#> 1 row(s) of entities have no grounding and can match nodes only by id.
#> Copied the values of entities to 4 node row(s), and added a protein-level row to 5 protein(s) with only PTM site rows.
protein_rows <- network$nodes[is.na(network$nodes$site), ]
knitr::kable(protein_rows[, c("entity_name", "node_role", "measured",
                              "logFC", "adj.pvalue")])
```

|     | entity_name | node_role      | measured |   logFC | adj.pvalue |
|:----|:------------|:---------------|:---------|--------:|-----------:|
| 1   | BRD4        | passed_cutoffs | TRUE     | -0.3303 |   0.026840 |
| 3   | LMNA        | passed_cutoffs | TRUE     |  0.4678 |   0.016000 |
| 5   | PML         | passed_cutoffs | TRUE     |  1.0720 |   0.000125 |
| 7   | IQGAP1      | passed_cutoffs | TRUE     |  0.7732 |   0.003143 |
| 9   | CDK18       | passed_cutoffs | TRUE     |  1.0550 |   0.019570 |
| 16  | MAPK3       | user_added     | TRUE     |  0.4066 |   0.308800 |
| 17  | AKT1        | user_added     | TRUE     |  0.4694 |   0.004076 |
| 18  | CDK1        | user_added     | FALSE    |      NA |         NA |
| 19  | CSNK2A1     | user_added     | FALSE    |      NA |         NA |
| 20  | MAPK1       | user_added     | TRUE     |  0.3307 |   0.532000 |
| 21  | SRC         | user_added     | FALSE    |      NA |         NA |
| 22  | PRKCA       | user_added     | FALSE    |      NA |         NA |
| 23  | PRKACA      | user_added     | TRUE     | -0.1899 |   0.151600 |

Kinases that were not measured in the global profiling experiment keep
`measured = FALSE`. Proteins that have sites in the network but no
protein-level row, such as HIF1A, were not measured in the global
profiling experiment either. Their sites are the `Adjusted = FALSE` rows
of `ADJUSTED.Model`, so their change is not adjusted for protein
abundance.

``` r

site_rows <- network$nodes[!is.na(network$nodes$site), ]
knitr::kable(site_rows[, c("entity_name", "site", "logFC", "adj.pvalue")])
```

|     | entity_name | site           |   logFC | adj.pvalue |
|:----|:------------|:---------------|--------:|-----------:|
| 2   | BRD4        | S1117          |  0.6929 |   0.034570 |
| 4   | LMNA        | S301           | -0.8042 |   0.041350 |
| 6   | PML         | S36_S38_S40    | -0.6853 |   0.042960 |
| 8   | IQGAP1      | Y1510          | -0.8199 |   0.034370 |
| 10  | CDK18       | S132           |  3.3970 |   0.001939 |
| 11  | PSMA3       | S250           |  0.6722 |   0.017440 |
| 12  | SHC1        | S139           |  0.7772 |   0.008324 |
| 13  | AQP2        | S256_S261_S264 | -1.5560 |   0.012470 |
| 14  | ADAM17      | T735           |  0.5357 |   0.013300 |
| 15  | HIF1A       | S451_T455      |  1.4060 |   0.011440 |

## Reading the network

With both levels on one network, each site can be read against its
protein. For example:

- PML is more abundant in tumor (protein `logFC` above 1), while its
  phosphorylation at S36/S38/S40, adjusted for that, is lower. MAPK1 and
  MAPK3 are known to phosphorylate these sites, and neither changed in
  abundance, so a change in their activity, not their amount, is one
  possible explanation.
- AKT1 is slightly more abundant in tumor, and its substrate site LMNA
  S301 is less phosphorylated after adjustment.

Edges are prior knowledge from INDRA, not evidence from this experiment.
Click an edge in the visualization below to read the evidence behind it.

## Network Visualization

`cytoscapeNetwork` draws a protein with sites as one node, with a small
node for each site around it. The protein node is colored by the
protein’s change in the global profiling experiment, and each site node
by the site’s adjusted change. The legend’s `logFC` covers both, so read
a protein’s color and its sites’ colors as changes from different
experiments. A protein with no global profiling result, such as HIF1A,
has no fill and a solid border. Kinases that were not measured have no
fill and a dashed border.

``` r

cytoscapeNetwork(network$nodes, network$edges,
                 displayLabelType = "entity_name")
```

To open the network in a browser instead, use `previewNetworkInBrowser`.

``` r

previewNetworkInBrowser(network$nodes, network$edges,
                        displayLabelType = "entity_name")
```

This package is distributed under the
[Artistic-2.0](https://opensource.org/licenses/Artistic-2.0) license.
However, its dependencies may have different licenses. In this example,
get_network with the INDRA backend depends on INDRA, which is
distributed under the [BSD
2-Clause](https://opensource.org/license/bsd-2-clause) license.
Furthermore, INDRA’s knowledge sources may have different licenses for
commercial applications. Please refer to the [INDRA
README](https://github.com/sorgerlab/indra?tab=readme-ov-file#indra-modules)
for more information on its knowledge sources and their associated
licenses.

## Session info

``` r

sessionInfo()
#> R version 4.6.1 (2026-06-24)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Ubuntu 24.04.5 LTS
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
#> LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
#> 
#> locale:
#>  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
#>  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
#>  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
#> [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
#> 
#> time zone: UTC
#> tzcode source: system (glibc)
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#> [1] MSstatsBioNet_1.5.5 MSstats_4.20.0      BiocStyle_2.40.0   
#> 
#> loaded via a namespace (and not attached):
#>  [1] tidyselect_1.2.1      viridisLite_0.4.3     dplyr_1.2.1          
#>  [4] farver_2.1.2          S7_0.2.2              bitops_1.1-0         
#>  [7] fastmap_1.2.0         XML_3.99-0.25         digest_0.6.39        
#> [10] lifecycle_1.0.5       survival_3.8-6        statmod_1.5.2        
#> [13] magrittr_2.0.5        compiler_4.6.1        r2r_0.1.2            
#> [16] rlang_1.3.0           sass_0.4.10           tools_4.6.1          
#> [19] yaml_2.3.12           data.table_1.18.6.1   knitr_1.52           
#> [22] stopwords_2.3         htmlwidgets_1.6.4     MSstatsConvert_1.22.1
#> [25] marray_1.90.0         xml2_1.6.0            RColorBrewer_1.1-3   
#> [28] KernSmooth_2.23-26    purrr_1.2.2           desc_1.4.3           
#> [31] grid_4.6.1            preprocessCore_1.74.0 caTools_1.18.4       
#> [34] ggplot2_4.0.3         scales_1.4.0          gtools_3.9.5         
#> [37] MASS_7.3-65           cli_3.6.6             crayon_1.5.3         
#> [40] rmarkdown_2.32        ragg_1.5.2            reformulas_0.4.4     
#> [43] generics_0.1.4        otel_0.2.0            httr_1.4.9           
#> [46] minqa_1.2.8           cachem_1.1.0          splines_4.6.1        
#> [49] parallel_4.6.1        BiocManager_1.30.27   vctrs_0.7.3          
#> [52] boot_1.3-32           Matrix_1.7-5          jsonlite_2.0.0       
#> [55] bookdown_0.48         ggrepel_0.9.8         systemfonts_1.3.2    
#> [58] limma_3.68.5          plotly_4.12.1         lgr_0.5.2            
#> [61] tidyr_1.3.2           jquerylib_0.1.4       glue_1.8.1           
#> [64] nloptr_2.2.1          pkgdown_2.2.1         gtable_0.3.6         
#> [67] lme4_2.0-6            mlapi_0.1.1           tibble_3.3.1         
#> [70] pillar_1.11.1         htmltools_0.5.9       gplots_3.3.0         
#> [73] float_0.3-4           rsparse_0.5.3         R6_2.6.1             
#> [76] textshaping_1.0.5     Rdpack_2.6.6          evaluate_1.0.5       
#> [79] lattice_0.22-9        rentrez_1.2.4         rbibutils_2.4.1      
#> [82] backports_1.5.1       RhpcBLASctl_0.23-42   bslib_0.12.0         
#> [85] text2vec_0.6.6        Rcpp_1.1.2            nlme_3.1-169         
#> [88] checkmate_2.3.4       xfun_0.62             fs_2.1.0             
#> [91] pkgconfig_2.0.3
```
