# PTM Analysis

## Installation

Run this code below to install MSstatsBioNet from bioconductor

``` r
if (!require("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

BiocManager::install("MSstatsBioNet")
```

### Dataset

We will be taking a subset of the dataset found in this
[paper](https://www.biorxiv.org/content/10.1101/2024.10.21.619348v1).
The table is the output of the MSstatsPTM function `groupComparisonPTM`
(filtered down to the columns that are actually needed)

``` r
input = data.table::fread(system.file(
    "extdata/garrido-2024.csv",
    package = "MSstatsBioNet"
))
head(input)
#>               Protein    Label     log2FC         SE    DF      pvalue
#>                <char>   <char>      <num>      <num> <int>       <num>
#> 1: P00533_S1039_S1042 t0 vs t1 -0.3200363 0.14994312     9 0.061580951
#> 2:       P00533_S1064 t0 vs t1  0.3566531 0.08915347     9 0.003108364
#> 3:   P00533_S991_S995 t0 vs t1 -0.1229037 0.10858298     9 0.286937655
#> 4:        P00533_T693 t0 vs t1 -0.0233444 0.17724459     9 0.898113182
#> 5:   P00533_T693_S695 t0 vs t1 -0.1659957 0.15000754     9 0.297173769
#> 6:       P00533_Y1110 t0 vs t1  0.2106324 0.09279031     9 0.049364434
#>    adj.pvalue  issue
#>         <num> <lgcl>
#> 1: 0.28024590     NA
#> 2: 0.06863598     NA
#> 3: 0.57907374     NA
#> 4: 0.96083634     NA
#> 5: 0.58809108     NA
#> 6: 0.25258914     NA
```

### ID Conversion

First, we need to convert the group comparison results to a format that
can be processed by INDRA. We can use the `annotateProteinInfoFromIndra`
function to obtain these mappings.

In the below example, we convert uniprot IDs to their corresponding Hgnc
IDs. We can also extract other information, such as hgnc gene name and
protein function.

``` r
library(MSstatsBioNet)
#> Loading required package: MSstats
#> 
#> Attaching package: 'MSstats'
#> The following object is masked from 'package:grDevices':
#> 
#>     savePlot
annotated_df = annotateProteinInfoFromIndra(input, "Uniprot")
head(annotated_df)
#>               Protein    Label     log2FC         SE    DF      pvalue
#>                <char>   <char>      <num>      <num> <int>       <num>
#> 1: P00533_S1039_S1042 t0 vs t1 -0.3200363 0.14994312     9 0.061580951
#> 2:       P00533_S1064 t0 vs t1  0.3566531 0.08915347     9 0.003108364
#> 3:   P00533_S991_S995 t0 vs t1 -0.1229037 0.10858298     9 0.286937655
#> 4:        P00533_T693 t0 vs t1 -0.0233444 0.17724459     9 0.898113182
#> 5:   P00533_T693_S695 t0 vs t1 -0.1659957 0.15000754     9 0.297173769
#> 6:       P00533_Y1110 t0 vs t1  0.2106324 0.09279031     9 0.049364434
#>    adj.pvalue  issue GlobalProtein UniprotId HgncId HgncName
#>         <num> <lgcl>        <char>    <char> <char>   <char>
#> 1: 0.28024590     NA        P00533    P00533   3236     EGFR
#> 2: 0.06863598     NA        P00533    P00533   3236     EGFR
#> 3: 0.57907374     NA        P00533    P00533   3236     EGFR
#> 4: 0.96083634     NA        P00533    P00533   3236     EGFR
#> 5: 0.58809108     NA        P00533    P00533   3236     EGFR
#> 6: 0.25258914     NA        P00533    P00533   3236     EGFR
#>    IsTranscriptionFactor IsKinase IsPhosphatase
#>                   <lgcl>   <lgcl>        <lgcl>
#> 1:                 FALSE     TRUE         FALSE
#> 2:                 FALSE     TRUE         FALSE
#> 3:                 FALSE     TRUE         FALSE
#> 4:                 FALSE     TRUE         FALSE
#> 5:                 FALSE     TRUE         FALSE
#> 6:                 FALSE     TRUE         FALSE
```

### Subnetwork Query

The package provides a function `getSubnetworkFromIndra` that retrieves
a subnetwork of proteins from the INDRA database based on differential
abundance analysis results. This function may help finding off target
subnetworks.

``` r
subnetwork <- getSubnetworkFromIndra(annotated_df, pvalueCutoff = 0.05, statement_types = c("Phosphorylation"), logfc_cutoff = 1, force_include_other = c("HGNC:3236"))
#> Warning in getSubnetworkFromIndra(annotated_df, pvalueCutoff = 0.05, statement_types = c("Phosphorylation"), : NOTICE: This function includes third-party software components
#>         that are licensed under the BSD 2-Clause License. Please ensure to
#>         include the third-party licensing agreements if redistributing this
#>         package or utilizing the results based on this package.
#>         See the LICENSE file for more details.
head(subnetwork$nodes)
#>        id hgncName        Site      logFC adj.pvalue
#>    <char>   <char>      <char>      <num>      <num>
#> 1: P00533     EGFR S1039_S1042 -0.3200363 0.28024590
#> 2: P00533     EGFR       S1064  0.3566531 0.06863598
#> 3: P00533     EGFR   S991_S995 -0.1229037 0.57907374
#> 4: P00533     EGFR        T693 -0.0233444 0.96083634
#> 5: P00533     EGFR   T693_S695 -0.1659957 0.58809108
#> 6: P00533     EGFR       Y1110  0.2106324 0.25258914
head(subnetwork$edges)
#>   source target site     interaction evidenceCount paperCount
#> 1 Q13480 P00533 <NA> Phosphorylation             2          1
#> 2 P00533 P28482 <NA> Phosphorylation             6          1
#> 3 P00533 Q13480 Y317 Phosphorylation             1          1
#> 4 P00533 Q13480 <NA> Phosphorylation            24          1
#> 5 P00533 Q13480 Y659 Phosphorylation            12          1
#> 6 P28482 Q13480 S454 Phosphorylation             3          1
#>                                                                                 evidenceLink
#> 1 https://db.indra.bio/statements/from_agents?subject=4066@HGNC&object=3236@HGNC&format=html
#> 2 https://db.indra.bio/statements/from_agents?subject=3236@HGNC&object=6871@HGNC&format=html
#> 3 https://db.indra.bio/statements/from_agents?subject=3236@HGNC&object=4066@HGNC&format=html
#> 4 https://db.indra.bio/statements/from_agents?subject=3236@HGNC&object=4066@HGNC&format=html
#> 5 https://db.indra.bio/statements/from_agents?subject=3236@HGNC&object=4066@HGNC&format=html
#> 6 https://db.indra.bio/statements/from_agents?subject=6871@HGNC&object=4066@HGNC&format=html
#>                                                        sourceCounts
#> 1                                        {"sparser": 1, "reach": 1}
#> 2                           {"sparser": 2, "reach": 3, "rlimsp": 1}
#> 3                                                     {"rlimsp": 1}
#> 4 {"reach": 14, "sparser": 5, "trips": 1, "bel_lc": 1, "rlimsp": 3}
#> 5              {"pc": 1, "pe": 1, "hprd": 8, "psp": 1, "signor": 1}
#> 6                                          {"hprd": 2, "signor": 1}
#>            stmt_hash
#> 1 -16511170354231919
#> 2   4046712837468004
#> 3  16300905780198239
#> 4 -27721131241182418
#> 5  28792790420132101
#> 6   7148668566491587
```

### Network Visualization

Visualize the subnetwork on your browser

``` r
previewNetworkInBrowser(subnetwork$nodes, subnetwork$edges, displayLabelType = "hgncName")
```

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
#> [1] MSstatsBioNet_1.3.4 MSstats_4.18.1      BiocStyle_2.38.0   
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
#> [25] curl_7.0.0            MSstatsConvert_1.20.0 marray_1.88.0        
#> [28] xml2_1.5.2            RColorBrewer_1.1-3    KernSmooth_2.23-26   
#> [31] purrr_1.2.1           desc_1.4.3            grid_4.5.2           
#> [34] preprocessCore_1.72.0 caTools_1.18.3        log4r_0.4.4          
#> [37] ggplot2_4.0.2         scales_1.4.0          gtools_3.9.5         
#> [40] MASS_7.3-65           cli_3.6.5             crayon_1.5.3         
#> [43] rmarkdown_2.30        ragg_1.5.1            reformulas_0.4.4     
#> [46] generics_0.1.4        otel_0.2.0            httr_1.4.8           
#> [49] minqa_1.2.8           cachem_1.1.0          splines_4.5.2        
#> [52] parallel_4.5.2        BiocManager_1.30.27   vctrs_0.7.1          
#> [55] boot_1.3-32           Matrix_1.7-4          jsonlite_2.0.0       
#> [58] bookdown_0.46         ggrepel_0.9.7         systemfonts_1.3.2    
#> [61] limma_3.66.0          plotly_4.12.0         lgr_0.5.2            
#> [64] jquerylib_0.1.4       tidyr_1.3.2           glue_1.8.0           
#> [67] nloptr_2.2.1          pkgdown_2.2.0         gtable_0.3.6         
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
