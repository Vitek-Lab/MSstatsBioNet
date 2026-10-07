# Metabolomics Analysis

## Installation

Run this code below to install MSstatsBioNet from bioconductor

``` r

if (!require("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

BiocManager::install("MSstatsBioNet")
```

## Purpose of this vignette

`MSstatsBioNet` is not limited to proteins. This vignette shows how to
take differential abundance results for small molecules (metabolites,
drugs, and drug metabolites), ground the compound names to database
identifiers, and place them in a network of prior knowledge from INDRA
together with the proteins they are known to interact with.

## Dataset

We will be taking a subset of the dataset found in this
[paper](https://doi.org/10.1371/journal.pone.0271794) (Panitchpakdi et
al., PLOS ONE 2022). Healthy volunteers took a single oral dose of the
antihistamine diphenhydramine, and plasma and skin swabs were sampled
over a time course and profiled by untargeted LC-MS/MS metabolomics. The
raw data are public on MassIVE under accession
[MSV000085944](https://massive.ucsd.edu/ProteoSAFe/dataset.jsp?accession=MSV000085944).

The table is the output of the MSstats function `groupComparison` for
the comparison of 480 minutes after the dose against baseline
(`480m vs 0m`). MSstats names the analyte column `Protein`. Since each
analyte here is a compound, the column has been renamed to `Analyte`.
The subset holds 25 compounds: diphenhydramine and its metabolites, a
few endogenous metabolites such as bile acids, and some unannotated
features.

``` r

input = data.table::fread(system.file(
    "extdata/panitchpakdi-2022.csv",
    package = "MSstatsBioNet"
))
knitr::kable(head(input[, .(Analyte, Label, log2FC, adj.pvalue, issue)]))
```

| Analyte | Label | log2FC | adj.pvalue | issue |
|:---|:---|---:|---:|:---|
| 1,7-bis(4-hydroxyphenyl)-5-\[(2R,3R,4S,5S,6R)-3,4,5-trihydroxy-6-(hydroxymethyl)oxan-2-yl\]oxyheptan-3-one | 480m vs 0m | Inf | 0.0000000 | oneConditionMissing |
| 19-Hydroxy-PGF | 480m vs 0m | -Inf | 0.0000000 | oneConditionMissing |
| 2-\[(2-carboxyethyl)amino\]propanoic acid | 480m vs 0m | -Inf | 0.0000000 | oneConditionMissing |
| 2-Hexaprenyl-6-Methoxy-3-Methyl-1,4-Benzoquinone | 480m vs 0m | Inf | 0.0000000 | oneConditionMissing |
| 253.6432_4.32 | 480m vs 0m | 1.998580 | 0.0045833 |  |
| 330.2274_2.75 | 480m vs 0m | 3.513129 | 0.0003918 |  |

Compounds named after the drug have some of the largest fold changes, as
expected 8 hours after the dose.

``` r

input[grepl("diphenhydramine|DIP_", Analyte, ignore.case = TRUE),
      .(Analyte, log2FC, adj.pvalue)]
#>                                                 Analyte   log2FC   adj.pvalue
#>                                                  <char>    <num>        <num>
#> 1: BT-HUMAN-Step2: Diphenhydramine [CARBONYL_REDUCTION] 7.528641 1.037161e-11
#> 2:                                  DIP_p_272.1644_16.5 4.521963 0.000000e+00
#> 3:                                      Diphenhydramine 5.419145 3.978151e-13
#> 4:                          N-desmethyl-diphenhydramine 6.454086 6.961764e-12
```

Rows with `issue == "oneConditionMissing"` were detected in only one of
the two time points, so their `log2FC` is `Inf` or `-Inf`.

``` r

table(input$issue, useNA = "ifany")
#> 
#>                     oneConditionMissing 
#>                  11                  14
```

## Building the entity table

The network functions work on an entity table with one row per analyte.
`prepare_entities` builds it from the `groupComparison` output. It needs
the column holding the analyte identifiers (`id_column`, `"Protein"` by
default, so here `"Analyte"`), the type of each analyte (`entity_type`),
and the identifier system of that column (`id_type`). Here the
identifiers are compound names, so `id_type = "chemical_name"`.

Both arguments take either one value for every row, or the name of a
column holding one value per row. We mark diphenhydramine as a `"drug"`
and every other compound as a `"metabolite"`, so that the drug gets its
own shape in the network.

``` r

library(MSstatsBioNet)
#> Loading required package: MSstats
#> 
#> Attaching package: 'MSstats'
#> The following object is masked from 'package:grDevices':
#> 
#>     savePlot
input$entity_type <- ifelse(input$Analyte == "Diphenhydramine",
                            "drug", "metabolite")
entities <- prepare_entities(input,
                             id_column = "Analyte",
                             entity_type = "entity_type",
                             id_type = "chemical_name")
knitr::kable(head(entities[, c("id", "entity_type", "id_type", "logFC",
                               "adj.pvalue")]))
```

| id | entity_type | id_type | logFC | adj.pvalue |
|:---|:---|:---|---:|---:|
| 1,7-bis(4-hydroxyphenyl)-5-\[(2R,3R,4S,5S,6R)-3,4,5-trihydroxy-6-(hydroxymethyl)oxan-2-yl\]oxyheptan-3-one | metabolite | chemical_name | Inf | 0.0000000 |
| 19-Hydroxy-PGF | metabolite | chemical_name | -Inf | 0.0000000 |
| 2-\[(2-carboxyethyl)amino\]propanoic acid | metabolite | chemical_name | -Inf | 0.0000000 |
| 2-Hexaprenyl-6-Methoxy-3-Methyl-1,4-Benzoquinone | metabolite | chemical_name | Inf | 0.0000000 |
| 253.6432_4.32 | metabolite | chemical_name | 1.998580 | 0.0045833 |
| 330.2274_2.75 | metabolite | chemical_name | 3.513129 | 0.0003918 |

The table keeps every compound, not only the significant ones, so that
the network can later recognize which of its nodes were measured.

## ID Conversion

INDRA identifies a compound by a database identifier such as a ChEBI ID,
not by its name. A backend converts names to identifiers.
[`indra_backend()`](https://vitek-lab.github.io/MSstatsBioNet/reference/indra_backend.md)
creates a backend for INDRA, and `backend_capabilities` lists the
conversions it supports.

``` r

indra <- indra_backend()
backend_capabilities(indra)$id_conversions
#> $protein
#> [1] "uniprot"          "uniprot_mnemonic" "hgnc_symbol"     
#> 
#> $ptm_site
#> [1] "uniprot"          "uniprot_mnemonic" "hgnc_symbol"     
#> 
#> $metabolite
#> [1] "chemical_name"
#> 
#> $lipid
#> [1] "chemical_name"
#> 
#> $drug
#> [1] "chemical_name"
```

For metabolites, lipids, and drugs, INDRA’s backend grounds each
`chemical_name` with [Gilda](https://github.com/gyorilab/gilda), INDRA’s
grounding service, and keeps whichever namespace Gilda returns (CHEBI,
MESH, PUBCHEM, …). `convert_ids` fills in the `namespace`, `entity_id`,
and `entity_name` columns.

So that this vignette builds without network access, the chunks that
call Gilda and INDRA (`convert_ids` and `get_network`) are not run when
the vignette is built. Their results were saved with
`inst/script/panitchpakdi-2022-network.R` and are loaded from
`extdata/panitchpakdi-2022-network.rds` instead. The output below comes
from those saved results.

``` r

entities <- convert_ids(indra, entities)
```

``` r

knitr::kable(entities[, c("id", "namespace", "entity_id", "entity_name")])
```

| id | namespace | entity_id | entity_name |
|:---|:---|:---|:---|
| 1,7-bis(4-hydroxyphenyl)-5-\[(2R,3R,4S,5S,6R)-3,4,5-trihydroxy-6-(hydroxymethyl)oxan-2-yl\]oxyheptan-3-one | CHEBI | CHEBI:183974 | 1,7-bis(4-hydroxyphenyl)-5-\[(2r,3r,4s,5s,6r)-3,4,5-trihydroxy-6-(hydroxymethyl)oxan-2-yl\]oxyheptan-3-one |
| 19-Hydroxy-PGF | MESH | C025183 | 19-hydroxyprostaglandin F |
| 2-\[(2-carboxyethyl)amino\]propanoic acid | CHEBI | CHEBI:15337 | beta-alanopine |
| 2-Hexaprenyl-6-Methoxy-3-Methyl-1,4-Benzoquinone | CHEBI | CHEBI:28711 | 2-hexaprenyl-6-methoxy-3-methyl-1,4-benzoquinone |
| 253.6432_4.32 | NA | NA | NA |
| 330.2274_2.75 | NA | NA | NA |
| alpha-Guaiaconic acid | CHEBI | CHEBI:193019 | alpha-guaiaconic acid |
| Batzelladine J | MESH | C505318 | batzelladine J |
| BT-HUMAN-Step2: Diphenhydramine \[CARBONYL_REDUCTION\] | NA | NA | NA |
| Dexpanthenol | CHEBI | CHEBI:27373 | pantothenol |
| DIP_p_272.1644_16.5 | NA | NA | NA |
| Diphenhydramine | CHEBI | CHEBI:4636 | diphenhydramine |
| Fluorene | CHEBI | CHEBI:28266 | fluorene |
| Galabiose | CHEBI | CHEBI:156273 | galabiose |
| GLYCOCHENODEOXYCHOLIC ACID | CHEBI;CHEBI | CHEBI:36274;CHEBI:180956 | glycochenodeoxycholic acid;Chenodeoxycholylglycine |
| Glycocholic acid | CHEBI | CHEBI:17687 | glycocholic acid |
| melibiose | CHEBI | CHEBI:28053 | melibiose |
| N-(2-Hydroxyethyl)octanamide | CHEBI | CHEBI:85302 | N-(octanoyl)ethanolamine |
| N-desmethyl-diphenhydramine | CHEBI | CHEBI:188299 | N-Desmethyldiphenhydramine |
| N-lactoyl-phenylalanine | MESH | C000723769 | N-lactoyl-phenylalanine |
| N-Lauroylethanolamine | CHEBI | CHEBI:85263 | N-(dodecanoyl)ethanolamine |
| Ryanodanol | MESH | C523743 | ryanodanol |
| talaromycin a | CHEBI | CHEBI:169379 | Talaromycin A |
| TAUROCHENODEOXYCHOLIC ACID | CHEBI | CHEBI:16525 | taurochenodeoxycholic acid |
| Trienediol | CHEBI | CHEBI:217780 | Trienediol |

A few things to note in the output:

- Diphenhydramine grounds to `CHEBI:4636`, and its N-demethylated
  metabolite grounds to its own ChEBI entry.
- Features named only by mass and retention time (e.g. `330.2274_2.75`)
  and predicted biotransformation products
  (e.g. `BT-HUMAN-Step2: Diphenhydramine [CARBONYL_REDUCTION]`) have no
  grounding (`NA`). They stay in the table but will be left out of the
  network.
- A name can ground to more than one candidate. Glycochenodeoxycholic
  acid matches two ChEBI entries, so its `namespace`, `entity_id`, and
  `entity_name` hold semicolon-joined values in the same order. The
  network query uses all candidates.

``` r

table(grounded = !is.na(entities$entity_id))
#> grounded
#> FALSE  TRUE 
#>     4    21
```

## Selecting Compounds

`select_entities` flags the compounds to send to INDRA by setting
`included_in_query`. Here we keep compounds with an adjusted p-value
below 0.05.

``` r

selected <- select_entities(entities, pvalue_cutoff = 0.05)
table(selected$included_in_query)
#> 
#> FALSE  TRUE 
#>    14    11
```

Compounds with an infinite fold change are never selected by the
cutoffs, because their p-values are not computed from a model. We come
back to them below.

## Network Query

### Measured compounds only

`get_network` sends the selected compounds to the backend and returns
`nodes` and `edges` tables.
[`subnetwork_query()`](https://vitek-lab.github.io/MSstatsBioNet/reference/subnetwork_query.md)
asks how the selected compounds are connected to each other, with no
other nodes added. It prints this question as a message.

``` r

network <- get_network(indra, selected, subnetwork_query())
#> Dropping 4 row(s) with no entity grounding (NA entity_id).
#> INDRA subnetwork: how are 7 selected entities connected to each other, with no other nodes added?
```

``` r

knitr::kable(network$nodes[, c("id", "entity_name", "namespace",
                               "entity_id", "logFC")])
```

| id | entity_name | namespace | entity_id | logFC |
|:---|:---|:---|:---|---:|
| GLYCOCHENODEOXYCHOLIC ACID | glycochenodeoxycholic acid;Chenodeoxycholylglycine | CHEBI;CHEBI | CHEBI:36274;CHEBI:180956 | 1.704549 |
| Glycocholic acid | glycocholic acid | CHEBI | CHEBI:17687 | 1.746646 |

``` r

knitr::kable(network$edges[, c("source", "target", "interaction",
                               "evidence_count")])
```

| source | target | interaction | evidence_count |
|:---|:---|:---|---:|
| Glycocholic acid | GLYCOCHENODEOXYCHOLIC ACID | Complex | 2 |
| GLYCOCHENODEOXYCHOLIC ACID | Glycocholic acid | Complex | 2 |

The literature rarely describes one small molecule acting on another
directly, so the network is small. Here only the two glycine-conjugated
bile acids are connected.

### Adding proteins with `include_entities`

Compounds are usually connected through proteins: the enzymes that make
or break them down, the transporters that move them, and the receptors
they act on. These proteins are not measured in a metabolomics
experiment, but they can be added to the query with the
`include_entities` argument of `get_network`, using
`"namespace:identifier"` IDs. Here we add:

- `HGNC:5182` (HRH1), the histamine H1 receptor that diphenhydramine
  blocks,
- `HGNC:2625` (CYP2D6), the main cytochrome P450 that metabolizes
  diphenhydramine,
- `HGNC:7967` (NR1H4, also known as FXR), the nuclear receptor for bile
  acids.

``` r

proteins <- c("HGNC:5182", "HGNC:2625", "HGNC:7967")
```

``` r

network <- get_network(indra, selected, subnetwork_query(),
                       include_entities = proteins)
#> Dropping 4 row(s) with no entity grounding (NA entity_id).
#> INDRA subnetwork: how are 7 selected entities and 3 added entities connected to each other, with no other nodes added?
```

``` r

knitr::kable(network$nodes[, c("id", "entity_type", "measured", "logFC")])
```

| id                         | entity_type | measured |    logFC |
|:---------------------------|:------------|:---------|---------:|
| Diphenhydramine            | drug        | TRUE     | 5.419145 |
| GLYCOCHENODEOXYCHOLIC ACID | metabolite  | TRUE     | 1.704549 |
| Glycocholic acid           | metabolite  | TRUE     | 1.746646 |
| TAUROCHENODEOXYCHOLIC ACID | metabolite  | TRUE     | 3.418690 |
| CYP2D6                     | protein     | FALSE    |       NA |
| NR1H4                      | protein     | FALSE    |       NA |
| HRH1                       | protein     | FALSE    |       NA |

``` r

knitr::kable(head(network$edges[, c("source", "target", "interaction",
                                    "evidence_count")], 10))
```

| source | target | interaction | evidence_count |
|:---|:---|:---|---:|
| NR1H4 | TAUROCHENODEOXYCHOLIC ACID | Activation | 1 |
| Glycocholic acid | NR1H4 | DecreaseAmount | 2 |
| Diphenhydramine | HRH1 | Inhibition | 41 |
| TAUROCHENODEOXYCHOLIC ACID | NR1H4 | Inhibition | 2 |
| TAUROCHENODEOXYCHOLIC ACID | NR1H4 | DecreaseAmount | 1 |
| NR1H4 | CYP2D6 | DecreaseAmount | 1 |
| NR1H4 | Glycocholic acid | Inhibition | 1 |
| Glycocholic acid | GLYCOCHENODEOXYCHOLIC ACID | Complex | 2 |
| Glycocholic acid | NR1H4 | Activation | 5 |
| NR1H4 | Glycocholic acid | Complex | 1 |

The added proteins have `measured = FALSE` and no statistics, because
they match no row of the entity table. The network now links
diphenhydramine to HRH1 and CYP2D6, and the bile acids to FXR.

`validate_network` checks that the network meets the contract that the
visualization and filtering functions expect. It returns the network
invisibly when it does, and gives an error listing the problems when it
does not.

``` r

validate_network(network)
```

### Compounds detected at one time point only

Set `include_infinite_fc = TRUE` in `select_entities` to also select
compounds with an infinite fold change. Note that a compound detected at
only one time point may be at the detection limit at the other, so treat
these compounds with more caution.

``` r

selected_inf <- select_entities(entities, pvalue_cutoff = 0.05,
                                include_infinite_fc = TRUE)
table(selected_inf$included_in_query)
#> 
#> TRUE 
#>   25
```

``` r

network_inf <- get_network(indra, selected_inf, subnetwork_query(),
                           include_entities = proteins)
#> Dropping 4 row(s) with no entity grounding (NA entity_id).
#> INDRA subnetwork: how are 21 selected entities and 3 added entities connected to each other, with no other nodes added?
```

``` r

setdiff(network_inf$nodes$id, network$nodes$id)
#> [1] "melibiose"
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

## Network Visualization

`cytoscapeNetwork` draws the network. The shape of a node comes from its
`entity_type`: diphenhydramine is a diamond, the other metabolites are
hexagons, and proteins are rounded rectangles. Proteins added through
`include_entities` have no fill and a dashed border, because they were
not measured. Click an edge to open its INDRA evidence.

``` r

cytoscapeNetwork(network$nodes, network$edges,
                 displayLabelType = "entity_name")
```

To open the network in a browser instead, use `previewNetworkInBrowser`.

``` r

previewNetworkInBrowser(network$nodes, network$edges,
                        displayLabelType = "entity_name")
```

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
#> [1] MSstatsBioNet_1.5.4 MSstats_4.20.0      BiocStyle_2.40.0   
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
#> [88] checkmate_2.3.4       xfun_0.61             fs_2.1.0             
#> [91] pkgconfig_2.0.3
```
