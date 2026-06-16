# Bootstrap the topic decomposition to find each topic's robust top words

Refits the NMF topic model on many bootstrap resamples of the papers and
reports, for every topic, how reliably each word stays among the topic's
top terms. This separates words that genuinely characterise a topic from
words that only surface in a single lucky fit, and lets you see how the
top-word lists change with and without the PPI view (run once per
`include_ppi`).

## Usage

``` r
bootstrapTopicModels(
  subnetwork,
  n_boot = 50,
  n_topics = 5,
  include_ppi = TRUE,
  n_top_terms = 10,
  min_term_count = 2,
  max_iter = 200,
  tol = 1e-04,
  seed = 1
)
```

## Arguments

- subnetwork:

  list with `nodes` and `edges` data.frames, e.g. the output of
  [`getSubnetworkFromIndra`](https://vitek-lab.github.io/MSstatsBioNet/reference/getSubnetworkFromIndra.md).

- n_boot:

  number of bootstrap resamples. Default 50.

- n_topics:

  number of topics (rank of the factorization). Default 5.

- include_ppi:

  logical; factorize the PPI/edge view jointly with the text (`TRUE`,
  default) or use paper words only (`FALSE`). See
  [`decomposeSubnetworkByTopic`](https://vitek-lab.github.io/MSstatsBioNet/reference/decomposeSubnetworkByTopic.md).

- n_top_terms:

  number of top words that define a topic's "top list" in each resample
  (the cutoff for the selection-frequency tally). Default 10.

- min_term_count:

  minimum corpus frequency for a word to be kept when building the text
  matrix. Default 2.

- max_iter:

  maximum number of NMF multiplicative-update iterations. Default 200.

- tol:

  relative-change tolerance for NMF early stopping. Default 1e-4.

- seed:

  random seed for the reference fit, the resampling, and each bootstrap
  NMF. Default 1.

## Value

A list with

- include_ppi, n_boot, n_topics:

  the settings used.

- topTerms:

  named list `topic_1` ... `topic_k`. Each is a data.frame sorted by
  `selection_freq`, with columns `term`, `selection_freq` (fraction of
  resamples the word was in this topic's top `n_top_terms`), and
  `mean_weight` (mean within-topic word weight across resamples). A word
  with `selection_freq` near 1 is a stable signature of the topic.

- reference:

  named list of the top `n_top_terms` words per topic from the single
  full-data fit, for comparison.

## Details

Because topic indices are arbitrary across fits (label switching), each
resample's topics are first aligned to a reference fit on the full data
by cosine similarity of their topic-word vectors. The papers are
resampled with replacement; the word vocabulary is held fixed (from the
full data) so topics remain comparable across resamples, and the NMF
seed is held fixed so the variability reported reflects *data*
resampling rather than random initialization.

## Note

**Beta feature:** This function is experimental and the API may change
without notice in future versions.

## See also

[`decomposeSubnetworkByTopic`](https://vitek-lab.github.io/MSstatsBioNet/reference/decomposeSubnetworkByTopic.md),
[`compareTopicModels`](https://vitek-lab.github.io/MSstatsBioNet/reference/compareTopicModels.md)

## Examples

``` r
if (FALSE) { # \dontrun{
input <- data.table::fread(system.file(
    "extdata/groupComparisonModel.csv",
    package = "MSstatsBioNet"
))
subnetwork <- getSubnetworkFromIndra(input)

# Top words with PPIs included vs. words only:
boot_ppi  <- bootstrapTopicModels(subnetwork, include_ppi = TRUE)
boot_text <- bootstrapTopicModels(subnetwork, include_ppi = FALSE)

head(boot_ppi$topTerms$topic_1)     # robust signature words for topic 1
boot_text$topTerms$topic_1
} # }
```
