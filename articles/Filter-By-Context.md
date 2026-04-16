# Filtering Subnetworks by Biological Context

## Overview

This vignette demonstrates how to use
[`filterSubnetworkByContext()`](https://vitek-lab.github.io/MSstatsBioNet/reference/filterSubnetworkByContext.md)
to filter a protein interaction subnetwork by the contextual relevance
of its supporting literature. The function:

1.  Retrieves evidence sentences from the INDRA database for each edge
    in the network
2.  Fetches the corresponding PubMed abstracts
3.  Scores each abstract against a user-supplied text query using TF-IDF
    cosine score
4.  Returns only the nodes, edges, and evidence whose abstracts exceed a
    score cutoff

This is useful when a subnetwork contains many edges supported by
literature from unrelated biological contexts, and you want to focus on
edges relevant to a specific research question — in this case, **DNA
damage repair in cancer**.

## Input Data

[`filterSubnetworkByContext()`](https://vitek-lab.github.io/MSstatsBioNet/reference/filterSubnetworkByContext.md)
expects a nodes and edges dataframe, typically produced by
[`getSubnetworkFromIndra()`](https://vitek-lab.github.io/MSstatsBioNet/reference/getSubnetworkFromIndra.md).
For this example we construct a small representative input table
directly, mimicking the structure of a proteomics experiment centred on
the DNA damage response kinase **CHK1**.

The input table contains one row per protein with columns for the
UniProt mnemonic identifier, the log2 fold-change, and the adjusted
p-value from a differential expression analysis.

``` r
input <- data.frame(
  Protein   = c("CHK1_HUMAN", "RFA1_HUMAN", "CLH1_HUMAN", "CRTC3_HUMAN"),
  log2FC    = c(2.31, 1.87, 1.45, 1.12),
  adj.pvalue = c(0.0021, 0.0089, 0.0310, 0.0490),
  stringsAsFactors = FALSE
)

input
```

    ##       Protein log2FC adj.pvalue
    ## 1  CHK1_HUMAN   2.31     0.0021
    ## 2  RFA1_HUMAN   1.87     0.0089
    ## 3  CLH1_HUMAN   1.45     0.0310
    ## 4 CRTC3_HUMAN   1.12     0.0490

All four proteins are up-regulated (positive log2FC) and statistically
significant (adj.pvalue \< 0.05).

------------------------------------------------------------------------

## Building the Subnetwork

### Step 1 — Annotate proteins with INDRA metadata

[`annotateProteinInfoFromIndra()`](https://vitek-lab.github.io/MSstatsBioNet/reference/annotateProteinInfoFromIndra.md)
maps UniProt mnemonics to HGNC gene identifiers and other metadata used
downstream by the INDRA query engine.

``` r
library(MSstatsBioNet)
annotated_df <- annotateProteinInfoFromIndra(input, "Uniprot_Mnemonic")
```

### Step 2 — Retrieve the interaction subnetwork

[`getSubnetworkFromIndra()`](https://vitek-lab.github.io/MSstatsBioNet/reference/getSubnetworkFromIndra.md)
queries the INDRA database for curated causal interactions among the
annotated proteins and returns a list containing `$nodes` and `$edges`
dataframes.

Key parameters used here:

- **`pvalueCutoff = 0.2`** — relaxed threshold to retain more candidate
  edges for downstream context filtering
- **`evidence_count_cutoff = 1`** — keep edges supported by at least one
  literature statement
- **`force_include_other = "HGNC:1925"`** — always include CHK1
  (HGNC:1925) regardless of significance, as it is the focal protein of
  interest
- **`filter_by_curation = FALSE`** — include both curated and
  automatically extracted interactions

``` r
subnetwork <- getSubnetworkFromIndra(
  annotated_df,
  pvalueCutoff          = 0.2,
  logfc_cutoff          = NULL,
  evidence_count_cutoff = 1,
  sources_filter        = NULL,
  force_include_other   = "HGNC:1925",
  filter_by_curation    = FALSE
)

# Inspect the unfiltered network
nrow(subnetwork$nodes)
nrow(subnetwork$edges)
```

------------------------------------------------------------------------

## Filtering by Context: Tag Count

### Defining the Query

The query string is compared against each PubMed abstract supporting the
network edges. A richer query — one that includes synonyms,
abbreviations, and related terms — improves recall under TF-IDF, which
relies on exact token matching rather than semantic understanding.

The expanded query below was produced with the help of a chatbot and
covers the major vocabulary used in the DNA damage repair and cancer
literature.

``` r
tags <- c(
  "dna damage repair",
  "cancer",
  "oncology",
  "dna repair",
  "genome integrity",
  "genomic instability",
  "double strand_break",
  "dsb",
  "single strand_break",
  "ssb",
  "base excision repair",
  "ber",
  "nucleotide excision repair",
  "ner",
  "mismatch repair",
  "mmr",
  "homologous recombination",
  "hr",
  "non homologous end joining",
  "nhej",
  "brca1",
  "brca2",
  "atm",
  "atr",
  "p53",
  "tp53",
  "parp",
  "tumor suppressor",
  "oncogene",
  "carcinogenesis",
  "tumorigenesis",
  "chemotherapy resistance",
  "radiation resistance",
  "genotoxic stress",
  "replication stress",
  "oxidative dna_damage",
  "somatic mutation",
  "tumor mutational burden",
  "tmb"
)
```

> **Tip:** You can iteratively refine `tags` by inspecting the scores in
> `filtered_network$evidence` and adding terms that appear frequently in
> high-scoring abstracts but are absent from your query.

[`filterSubnetworkByContext()`](https://vitek-lab.github.io/MSstatsBioNet/reference/filterSubnetworkByContext.md)
ties everything together. The `cutoff` parameter controls stringency —
only edges whose supporting abstracts score at or above this value are
retained.

``` r
filtered_network <- filterSubnetworkByContext(
  nodes             = subnetwork$nodes,
  edges             = subnetwork$edges,
  method            = "tag_count",
  cutoff = 3,
  query             = tags
)
```

The function prints a progress summary to the console:

    Processing N unique statement hashes...
    Fetching M abstracts...
    Progress: M/M (100.0%)
    Done fetching abstracts!

    X / M abstracts passed score cutoff (>= 0.10)
    Retained: A edges (of B), C nodes (of D), E evidence rows (of F)

### Filtered nodes

``` r
filtered_network$nodes
```

Only proteins connected by at least one contextually relevant edge are
retained.

### Filtered edges

``` r
filtered_network$edges
```

Each row represents a causal interaction (e.g. phosphorylation,
activation) supported by literature that passed the score threshold.

### Evidence with scores

``` r
filtered_network$evidence
```

The evidence dataframe contains the following columns:

| Column         | Description                                  |
|----------------|----------------------------------------------|
| `source`       | Source protein / gene                        |
| `target`       | Target protein / gene                        |
| `interaction`  | Interaction type (e.g. Phosphorylation)      |
| `site`         | Modification site if applicable              |
| `evidenceLink` | URL to the INDRA evidence viewer             |
| `stmt_hash`    | Unique INDRA statement identifier            |
| `text`         | Sentence extracted from the supporting paper |
| `pmid`         | PubMed ID of the source article              |
| `score`        | Cosine score of the abstract vs. query       |

You can sort by score to identify the most on-topic supporting evidence:

``` r
filtered_network$evidence[
  order(filtered_network$evidence$score, decreasing = TRUE), 
]
```

### Defining a cutoff

``` r
# Run with permissive cutoff to see full score distribution
exploratory <- filterSubnetworkByContext(
  nodes             = subnetwork$nodes,
  edges             = subnetwork$edges,
  cutoff = 0.0,
  query             = tags
)

summary(exploratory$evidence$score)
hist(exploratory$evidence$score,
     breaks = 30,
     main   = "Distribution of abstract scores",
     xlab   = "Number of tags matched",
     col    = "steelblue")
```

------------------------------------------------------------------------

## Filtering by Context: Cosine score

### Defining the Query

The query string is compared against each PubMed abstract supporting the
network edges. A richer query — one that includes synonyms,
abbreviations, and related terms — improves recall under TF-IDF, which
relies on exact token matching rather than semantic understanding.

The expanded query below was produced with the help of a chatbot and
covers the major vocabulary used in the DNA damage repair and cancer
literature.

``` r
my_query <- "DNA damage repair cancer oncology DNA repair genome integrity
  genomic instability double strand break DSB single strand break SSB
  base excision repair BER nucleotide excision repair NER mismatch repair MMR
  homologous recombination HR non-homologous end joining NHEJ BRCA1 BRCA2
  ATM ATR p53 TP53 PARP tumor suppressor oncogene carcinogenesis tumorigenesis
  chemotherapy resistance radiation resistance genotoxic stress replication stress
  oxidative DNA damage somatic mutation tumor mutational burden TMB"
```

> **Tip:** You can iteratively refine `my_query` by inspecting the
> scores in `filtered_network$evidence` and adding terms that appear
> frequently in high-scoring abstracts but are absent from your query.

[`filterSubnetworkByContext()`](https://vitek-lab.github.io/MSstatsBioNet/reference/filterSubnetworkByContext.md)
ties everything together. The `cutoff` parameter controls stringency —
only edges whose supporting abstracts score at or above this value are
retained.

``` r
filtered_network <- filterSubnetworkByContext(
  nodes             = subnetwork$nodes,
  edges             = subnetwork$edges,
  method            = "cosine",
  cutoff = 0.10,
  query             = my_query
)
```

The function prints a progress summary to the console:

    Processing N unique statement hashes...
    Fetching M abstracts...
    Progress: M/M (100.0%)
    Done fetching abstracts!

    X / M abstracts passed score cutoff (>= 0.10)
    Retained: A edges (of B), C nodes (of D), E evidence rows (of F)

### Filtered nodes

``` r
filtered_network$nodes
```

Only proteins connected by at least one contextually relevant edge are
retained.

### Filtered edges

``` r
filtered_network$edges
```

Each row represents a causal interaction (e.g. phosphorylation,
activation) supported by literature that passed the score threshold.

### Evidence with scores

``` r
filtered_network$evidence
```

The evidence dataframe contains the following columns:

| Column         | Description                                      |
|----------------|--------------------------------------------------|
| `source`       | Source protein / gene                            |
| `target`       | Target protein / gene                            |
| `interaction`  | Interaction type (e.g. Phosphorylation)          |
| `site`         | Modification site if applicable                  |
| `evidenceLink` | URL to the INDRA evidence viewer                 |
| `stmt_hash`    | Unique INDRA statement identifier                |
| `text`         | Sentence extracted from the supporting paper     |
| `pmid`         | PubMed ID of the source article                  |
| `score`        | Relevance score (tag count or cosine similarity) |

You can sort by score to identify the most on-topic supporting evidence:

``` r
filtered_network$evidence[
  order(filtered_network$evidence$score, decreasing = TRUE), 
]
```

### Choosing a score Cutoff

The right cutoff depends on how broadly the query overlaps with the
literature in your network. As a rough guide:

| Cutoff   | Effect                                                   |
|----------|----------------------------------------------------------|
| `0.05`   | Permissive — removes only completely off-topic abstracts |
| `0.10`   | Recommended default for domain-specific queries          |
| `0.20`   | Stringent — retains only highly on-topic edges           |
| `> 0.30` | Very stringent — use only with highly specific queries   |

To explore the score distribution before committing to a cutoff, run the
function at a low threshold and inspect the scores:

``` r
# Run with permissive cutoff to see full score distribution
exploratory <- filterSubnetworkByContext(
  nodes             = subnetwork$nodes,
  edges             = subnetwork$edges,
  cutoff = 0.0,
  method = "cosine",
  query             = my_query
)

summary(exploratory$evidence$score)
hist(exploratory$evidence$score,
     breaks = 30,
     main   = "Distribution of abstract scores",
     xlab   = "Cosine score to query",
     col    = "steelblue")
```

------------------------------------------------------------------------

## Session Info

``` r
sessionInfo()
#> R version 4.5.3 (2026-03-11)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Ubuntu 24.04.4 LTS
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
#> loaded via a namespace (and not attached):
#>  [1] digest_0.6.39     desc_1.4.3        R6_2.6.1          fastmap_1.2.0    
#>  [5] xfun_0.57         cachem_1.1.0      knitr_1.51        htmltools_0.5.9  
#>  [9] rmarkdown_2.31    lifecycle_1.0.5   cli_3.6.6         sass_0.4.10      
#> [13] pkgdown_2.2.0     textshaping_1.0.5 jquerylib_0.1.4   systemfonts_1.3.2
#> [17] compiler_4.5.3    tools_4.5.3       ragg_1.5.2        bslib_0.10.0     
#> [21] evaluate_1.0.5    yaml_2.3.12       otel_0.2.0        jsonlite_2.0.0   
#> [25] rlang_1.2.0       fs_2.0.1          htmlwidgets_1.6.4
```
