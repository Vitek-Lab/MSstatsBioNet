# MSstatsBioNet Introduction

## Installation

Run this code below to install MSstatsBioNet from bioconductor

``` r
if (!require("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

BiocManager::install("MSstatsBioNet")
```

## Purpose of MSstatsBioNet

The `MSstatsBioNet` package is a member of the `MSstats` family of
packages. It contains a set of functions for interpretation of mass
spectrometry (MS) statistical analysis results in the context of
protein-protein interaction networks. The package is designed to be used
in conjunction with the `MSstats` package.

## Dataset

We will be taking a subset of the dataset found in this
[paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC7331093/).

``` r
input = data.table::fread(system.file(
    "extdata/msstats.csv",
    package = "MSstatsBioNet"
))
```

## MSstats Convert from Upstream Dataset

``` r
library(MSstatsConvert)
msstats_imported = FragPipetoMSstatsFormat(input, use_log_file = FALSE)
#> INFO  [2026-03-11 18:26:05] ** Raw data from FragPipe imported successfully.
#> INFO  [2026-03-11 18:26:05] ** Using annotation extracted from quantification data.
#> INFO  [2026-03-11 18:26:05] ** Run labels were standardized to remove symbols such as '.' or '%'.
#> INFO  [2026-03-11 18:26:05] ** The following options are used:
#>   - Features will be defined by the columns: PeptideSequence, PrecursorCharge, FragmentIon, ProductCharge
#>   - Shared peptides will be removed.
#>   - Proteins with single feature will not be removed.
#>   - Features with less than 3 measurements across runs will be removed.
#> INFO  [2026-03-11 18:26:05] ** Features with all missing measurements across runs are removed.
#> INFO  [2026-03-11 18:26:05] ** Shared peptides are removed.
#> INFO  [2026-03-11 18:26:05] ** Multiple measurements in a feature and a run are summarized by summaryforMultipleRows: max
#> INFO  [2026-03-11 18:26:05] ** Features with one or two measurements across runs are removed.
#> INFO  [2026-03-11 18:26:05] ** Run annotation merged with quantification data.
#> INFO  [2026-03-11 18:26:05] ** Features with one or two measurements across runs are removed.
#> INFO  [2026-03-11 18:26:05] ** Fractionation handled.
#> INFO  [2026-03-11 18:26:05] ** Updated quantification data to make balanced design. Missing values are marked by NA
#> INFO  [2026-03-11 18:26:05] ** Finished preprocessing. The dataset is ready to be processed by the dataProcess function.
head(msstats_imported)
#>   ProteinName PeptideSequence PrecursorCharge FragmentIon ProductCharge
#> 1      P05023   AVAGDASESALLK               2         b12             1
#> 2      P05023   AVAGDASESALLK               2         b12             1
#> 3      P05023   AVAGDASESALLK               2         b12             1
#> 4      P05023   AVAGDASESALLK               2         b12             1
#> 5      P05023   AVAGDASESALLK               2         b12             1
#> 6      P05023   AVAGDASESALLK               2         b12             1
#>   IsotopeLabelType Condition BioReplicate
#> 1                L       NAT            1
#> 2                L         T            1
#> 3                L       NAT            2
#> 4                L         T            2
#> 5                L       NAT            3
#> 6                L         T            3
#>                                              Run Fraction Intensity
#> 1 CPTAC_CCRCC_W_JHU_20190112_LUMOS_C3L-00418_NAT        1  69915.91
#> 2   CPTAC_CCRCC_W_JHU_20190112_LUMOS_C3L-00418_T        1 101726.53
#> 3 CPTAC_CCRCC_W_JHU_20190112_LUMOS_C3L-00606_NAT        1  36215.18
#> 4   CPTAC_CCRCC_W_JHU_20190112_LUMOS_C3L-00606_T        1  24927.36
#> 5 CPTAC_CCRCC_W_JHU_20190112_LUMOS_C3L-01882_NAT        1        NA
#> 6   CPTAC_CCRCC_W_JHU_20190112_LUMOS_C3L-01882_T        1        NA
```

We will first convert the input data to a format that can be processed
by MSstats. In this example, the `MetamorpheusToMSstatsFormat` function
is used to convert the input data from Metamorpheus to MSstats format.
The MSstats format will contain the necessary information for downstream
analysis, such as peptide information, abundance values, run ID, and
experimental annotation information.

## MSstats Process and GroupComparison

``` r
library(MSstats)
#> 
#> Attaching package: 'MSstats'
#> The following object is masked from 'package:grDevices':
#> 
#>     savePlot
QuantData <- dataProcess(msstats_imported, use_log_file = FALSE)
#> INFO  [2026-03-11 18:26:07] ** Log2 intensities under cutoff = 9.2885  were considered as censored missing values.
#> INFO  [2026-03-11 18:26:07] ** Log2 intensities = NA were considered as censored missing values.
#> INFO  [2026-03-11 18:26:07] ** Use top100 features that have highest average of log2(intensity) across runs.
#> INFO  [2026-03-11 18:26:07] 
#>  # proteins: 10
#>  # peptides per protein: 1-8
#>  # features per peptide: 7-12
#> INFO  [2026-03-11 18:26:07] 
#>                     NAT T
#>              # runs   5 5
#>     # bioreplicates   5 5
#>  # tech. replicates   1 1
#> INFO  [2026-03-11 18:26:07] Some features are completely missing in at least one condition:  
#>  ELEAEIQQLR_2_b5_1,
#>  ELEAEIQQLR_2_b8_1,
#>  NLEAVETLGSTSTIC(UniMod:4)SDK_3_b13_2,
#>  NLEAVETLGSTSTIC(UniMod:4)SDK_3_b3_1,
#>  NLEAVETLGSTSTIC(UniMod:4)SDK_3_b4_1 ...
#> INFO  [2026-03-11 18:26:07]  == Start the summarization per subplot...
#>   |                                                                              |                                                                      |   0%  |                                                                              |=======                                                               |  10%  |                                                                              |==============                                                        |  20%  |                                                                              |=====================                                                 |  30%  |                                                                              |============================                                          |  40%  |                                                                              |===================================                                   |  50%  |                                                                              |==========================================                            |  60%  |                                                                              |=================================================                     |  70%  |                                                                              |========================================================              |  80%  |                                                                              |===============================================================       |  90%  |                                                                              |======================================================================| 100%
#> INFO  [2026-03-11 18:26:07]  == Summarization is done.
model <- groupComparison(
    contrast.matrix = "pairwise",
    data = QuantData,
    use_log_file = FALSE
)
#> INFO  [2026-03-11 18:26:08]  == Start to test and get inference in whole plot ...
#>   |                                                                              |                                                                      |   0%  |                                                                              |=======                                                               |  10%  |                                                                              |==============                                                        |  20%  |                                                                              |=====================                                                 |  30%  |                                                                              |============================                                          |  40%  |                                                                              |===================================                                   |  50%  |                                                                              |==========================================                            |  60%  |                                                                              |=================================================                     |  70%  |                                                                              |========================================================              |  80%  |                                                                              |===============================================================       |  90%  |                                                                              |======================================================================| 100%
#> INFO  [2026-03-11 18:26:08]  == Comparisons for all proteins are done.
head(model$ComparisonResult)
#>   Protein    Label     log2FC        SE     Tvalue DF       pvalue  adj.pvalue
#> 1  O00217 NAT vs T  2.0285031 0.4364177   4.648077  4 0.0096753522 0.013821932
#> 2  O00330 NAT vs T  1.3000941 0.1320659   9.844285  3 0.0022285170 0.004457034
#> 3  O60313 NAT vs T  0.9299641 0.2387992   3.894336  4 0.0176257617 0.019584180
#> 4  O60879 NAT vs T -1.9484511 0.1695323 -11.493098  4 0.0003271868 0.001635934
#> 5  O75306 NAT vs T  2.4745040 0.3530857   7.008225  4 0.0021825029 0.004457034
#> 6  P05023 NAT vs T  1.8391155 0.2122058   8.666660  4 0.0009753219 0.003251073
#>   issue MissingPercentage ImputationPercentage
#> 1    NA        0.14166667           0.14166667
#> 2    NA        0.20833333           0.10833333
#> 3    NA        0.29756098           0.29756098
#> 4    NA        0.09583333           0.09583333
#> 5    NA        0.16111111           0.16111111
#> 6    NA        0.11702128           0.11702128
```

Next, we will preprocess the data using the `dataProcess` function and
perform a statistical analysis using the `groupComparison` function. The
output of `groupComparison` is a table containing a list of
differentially abundant proteins with their associated p-values and log
fold changes.

## MSstatsBioNet Analysis

### ID Conversion

First, we need to convert the group comparison results to a format that
can be processed by INDRA. The `getSubnetworkFromIndra` function
requires a table containing HGNC IDs. We can use the
`annotateProteinInfoFromIndra` function to obtain these mappings.

In the below example, we convert uniprot IDs to their corresponding Hgnc
IDs. We can also extract other information, such as hgnc gene name and
protein function.

``` r
library(MSstatsBioNet)
annotated_df = annotateProteinInfoFromIndra(model$ComparisonResult, "Uniprot")
head(annotated_df)
#>   Protein    Label     log2FC        SE     Tvalue DF       pvalue  adj.pvalue
#> 1  O00217 NAT vs T  2.0285031 0.4364177   4.648077  4 0.0096753522 0.013821932
#> 2  O00330 NAT vs T  1.3000941 0.1320659   9.844285  3 0.0022285170 0.004457034
#> 3  O60313 NAT vs T  0.9299641 0.2387992   3.894336  4 0.0176257617 0.019584180
#> 4  O60879 NAT vs T -1.9484511 0.1695323 -11.493098  4 0.0003271868 0.001635934
#> 5  O75306 NAT vs T  2.4745040 0.3530857   7.008225  4 0.0021825029 0.004457034
#> 6  P05023 NAT vs T  1.8391155 0.2122058   8.666660  4 0.0009753219 0.003251073
#>   issue MissingPercentage ImputationPercentage GlobalProtein UniprotId HgncId
#> 1    NA        0.14166667           0.14166667        O00217    O00217   7715
#> 2    NA        0.20833333           0.10833333        O00330    O00330  21350
#> 3    NA        0.29756098           0.29756098        O60313    O60313   8140
#> 4    NA        0.09583333           0.09583333        O60879    O60879   2877
#> 5    NA        0.16111111           0.16111111        O75306    O75306   7708
#> 6    NA        0.11702128           0.11702128        P05023    P05023    799
#>   HgncName IsTranscriptionFactor IsKinase IsPhosphatase
#> 1   NDUFS8                 FALSE    FALSE         FALSE
#> 2     PDHX                 FALSE    FALSE         FALSE
#> 3     OPA1                 FALSE    FALSE         FALSE
#> 4   DIAPH2                 FALSE    FALSE         FALSE
#> 5   NDUFS2                 FALSE    FALSE         FALSE
#> 6   ATP1A1                 FALSE    FALSE         FALSE
```

### Subnetwork Query

The package provides a function `getSubnetworkFromIndra` that retrieves
a subnetwork of proteins from the INDRA database based on differential
abundance analysis results.

``` r
subnetwork <- getSubnetworkFromIndra(
    annotated_df, 
    pvalueCutoff = 0.05, 
    statement_types = c("Complex", "IncreaseAmount", "DecreaseAmount", "Inhibition", "Activation", "Phosphorylation")
)
#> Warning in getSubnetworkFromIndra(annotated_df, pvalueCutoff = 0.05, statement_types = c("Complex", : NOTICE: This function includes third-party software components
#>         that are licensed under the BSD 2-Clause License. Please ensure to
#>         include the third-party licensing agreements if redistributing this
#>         package or utilizing the results based on this package.
#>         See the LICENSE file for more details.
head(subnetwork$nodes)
#>       id hgncName Site     logFC  adj.pvalue
#> 1 O00217   NDUFS8 <NA> 2.0285031 0.013821932
#> 3 O60313     OPA1 <NA> 0.9299641 0.019584180
#> 5 O75306   NDUFS2 <NA> 2.4745040 0.004457034
#> 6 P05023   ATP1A1 <NA> 1.8391155 0.003251073
#> 7 P05067      APP <NA> 0.7360012 0.020306662
#> 8 P05090     APOD <NA> 0.5683951 0.013715050
head(subnetwork$edges)
#>   source target site interaction evidenceCount paperCount
#> 1 P05023 O75306 <NA>     Complex             1          1
#> 2 O75306 P08574 <NA>     Complex             1          1
#> 3 P05067 O60313 <NA>  Activation             2          1
#> 4 O60313 O00217 <NA>     Complex             1          1
#> 5 P05362 P05067 <NA>     Complex            13          1
#> 6 O75306 P05067 <NA>     Complex             1          1
#>                                                                                 evidenceLink
#> 1  https://db.indra.bio/statements/from_agents?subject=799@HGNC&object=7708@HGNC&format=html
#> 2 https://db.indra.bio/statements/from_agents?subject=7708@HGNC&object=2579@HGNC&format=html
#> 3  https://db.indra.bio/statements/from_agents?subject=620@HGNC&object=8140@HGNC&format=html
#> 4 https://db.indra.bio/statements/from_agents?subject=8140@HGNC&object=7715@HGNC&format=html
#> 5  https://db.indra.bio/statements/from_agents?subject=5344@HGNC&object=620@HGNC&format=html
#> 6  https://db.indra.bio/statements/from_agents?subject=7708@HGNC&object=620@HGNC&format=html
#>                  sourceCounts          stmt_hash
#> 1              {"biogrid": 1}  -5813063534036006
#> 2              {"biogrid": 1}   6349003830434161
#> 3                {"reach": 2}   3948742039105656
#> 4              {"biogrid": 1} -19747883270157675
#> 5 {"sparser": 10, "reach": 3} -20220236678417803
#> 6              {"biogrid": 1}  22463147519060585
```

This package is distributed under the
[Artistic-2.0](https://opensource.org/licenses/Artistic-2.0) license.
However, its dependencies may have different licenses. In this example,
getSubnetworkFromIndra depends on INDRA, which is distributed under the
[BSD 2-Clause](https://opensource.org/license/bsd-2-clause) license.
Furthermore, INDRA’s knowledge sources may have different licenses for
commercial applications. Please refer to the [INDRA
README](https://github.com/sorgerlab/indra?tab=readme-ov-file#indra-modules)
for more information on its knowledge sources and their associated
licenses.

### Visualize Networks

The function `previewNetworkInBrowser` then takes the output of
`getSubnetworkFromIndra` and visualizes the subnetwork. The function
requires an internet browser to view the subnetwork

``` r
previewNetworkInBrowser(subnetwork$nodes, subnetwork$edges, displayLabelType = "hgncName")
```

In the network diagram displayed using CytoscapeJS, you should see three
arrows connecting two nodes, APP and APOD. These arrows represent the
interactions between these two proteins, notably Activation,
IncreaseAmount (e.g. via transcriptional regulation), and forming a
protein complex.

## Session info

``` r
sessionInfo()
#> R version 4.5.2 (2025-10-31)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Ubuntu 24.04.3 LTS
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
#> [1] MSstatsBioNet_1.3.4   MSstats_4.18.1        MSstatsConvert_1.20.0
#> [4] BiocStyle_2.38.0     
#> 
#> loaded via a namespace (and not attached):
#>  [1] tidyselect_1.2.1      viridisLite_0.4.3     dplyr_1.2.0          
#>  [4] farver_2.1.2          S7_0.2.1              bitops_1.0-9         
#>  [7] fastmap_1.2.0         lazyeval_0.2.2        XML_3.99-0.22        
#> [10] digest_0.6.39         lifecycle_1.0.5       survival_3.8-3       
#> [13] statmod_1.5.1         magrittr_2.0.4        compiler_4.5.2       
#> [16] r2r_0.1.2             rlang_1.1.7           sass_0.4.10          
#> [19] tools_4.5.2           yaml_2.3.12           data.table_1.18.2.1  
#> [22] knitr_1.51            stopwords_2.3         htmlwidgets_1.6.4    
#> [25] curl_7.0.0            marray_1.88.0         xml2_1.5.2           
#> [28] RColorBrewer_1.1-3    KernSmooth_2.23-26    purrr_1.2.1          
#> [31] desc_1.4.3            grid_4.5.2            preprocessCore_1.72.0
#> [34] caTools_1.18.3        log4r_0.4.4           ggplot2_4.0.2        
#> [37] scales_1.4.0          gtools_3.9.5          MASS_7.3-65          
#> [40] cli_3.6.5             crayon_1.5.3          rmarkdown_2.30       
#> [43] ragg_1.5.1            reformulas_0.4.4      generics_0.1.4       
#> [46] otel_0.2.0            httr_1.4.8            minqa_1.2.8          
#> [49] cachem_1.1.0          splines_4.5.2         parallel_4.5.2       
#> [52] BiocManager_1.30.27   vctrs_0.7.1           boot_1.3-32          
#> [55] Matrix_1.7-4          jsonlite_2.0.0        bookdown_0.46        
#> [58] ggrepel_0.9.7         systemfonts_1.3.2     limma_3.66.0         
#> [61] plotly_4.12.0         lgr_0.5.2             jquerylib_0.1.4      
#> [64] tidyr_1.3.2           glue_1.8.0            nloptr_2.2.1         
#> [67] pkgdown_2.2.0         stringi_1.8.7         gtable_0.3.6         
#> [70] lme4_2.0-1            mlapi_0.1.1           tibble_3.3.1         
#> [73] pillar_1.11.1         htmltools_0.5.9       gplots_3.3.0         
#> [76] float_0.3-3           rsparse_0.5.3         R6_2.6.1             
#> [79] textshaping_1.0.5     Rdpack_2.6.6          evaluate_1.0.5       
#> [82] lattice_0.22-7        rentrez_1.2.4         rbibutils_2.4.1      
#> [85] backports_1.5.0       RhpcBLASctl_0.23-42   bslib_0.10.0         
#> [88] text2vec_0.6.6        Rcpp_1.1.1            nlme_3.1-168         
#> [91] checkmate_2.3.4       xfun_0.56             fs_1.6.7             
#> [94] pkgconfig_2.0.3
```
