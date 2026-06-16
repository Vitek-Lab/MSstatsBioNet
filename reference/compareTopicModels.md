# Test whether including PPIs changes topic structure beyond random chance

Quantifies how much the topic decomposition produced by
[`decomposeSubnetworkByTopic`](https://vitek-lab.github.io/MSstatsBioNet/reference/decomposeSubnetworkByTopic.md)
changes when the PPI/edge view is included (`include_ppi = TRUE`) versus
excluded (`include_ppi = FALSE`), and separates that change from the
run-to-run variability that NMF produces just from its random
initialization.

## Usage

``` r
compareTopicModels(
  subnetwork,
  seeds = seq_len(20),
  n_topics = 5,
  unit = c("edges", "papers"),
  min_term_count = 2,
  max_iter = 200,
  tol = 1e-04
)
```

## Arguments

- subnetwork:

  list with `nodes` and `edges` data.frames, e.g. the output of
  [`getSubnetworkFromIndra`](https://vitek-lab.github.io/MSstatsBioNet/reference/getSubnetworkFromIndra.md).

- seeds:

  integer vector of NMF seeds to fit (at least 2). Default `1:20`.

- n_topics:

  number of topics (rank of the factorization). Default 5.

- unit:

  either `"edges"` (compare edge-to-topic assignments, the default,
  matching the subnetworks the decomposition returns) or `"papers"`
  (compare paper-to-topic assignments).

- min_term_count:

  minimum corpus frequency for a word to be kept when building the text
  matrix. Default 2.

- max_iter:

  maximum number of NMF multiplicative-update iterations. Default 200.

- tol:

  relative-change tolerance for NMF early stopping. Default 1e-4.

## Value

A list with

- unit:

  the comparison unit used.

- seeds:

  the seeds fitted.

- n_topics:

  the effective number of topics.

- ari:

  list of numeric vectors `within_joint`, `within_text`, and `between`
  (matched seeds).

- summary:

  data.frame of median/mean ARI and count per comparison.

- test:

  the [`wilcox.test`](https://rdrr.io/r/stats/wilcox.test.html) object
  comparing the between distribution against the pooled within
  distributions (`alternative = "less"`).

- consensus:

  list of consensus (co-membership) matrices, `joint` and `text`, across
  seeds.

- dispersion:

  named numeric vector of consensus dispersion coefficients (1 =
  identical clustering across all seeds).

- partitions:

  list of the raw per-seed partitions, `joint` and `text`, for further
  inspection.

## Details

NMF converges to a local optimum that depends on the random seed, so a
single joint-vs-text comparison conflates the real effect of the PPI
view with optimization noise. This function instead refits both modes
across many seeds and compares three distributions of partition
agreement (Adjusted Rand Index, ARI):

- within_joint:

  ARI between pairs of joint runs (different seeds) — how much the joint
  solution wobbles on its own.

- within_text:

  ARI between pairs of text-only runs — the same for the text-only
  solution.

- between:

  ARI between the joint and text-only run at the *same* seed. Because
  both modes draw `W` and `H_text` from the same seeded stream, a
  matched seed gives both modes an identical initialization, so this
  isolates the effect of adding the PPI view from the starting point.

If the between-mode ARI is systematically lower than the within-mode
ARIs, the PPI view changes the topic structure more than chance would —
a one-sided Wilcoxon rank-sum test (`between < within`) puts a p-value
on it. If the between distribution sits inside the within distributions,
the apparent difference is just optimization noise.

The expensive, network-bound steps (evidence extraction, abstract
fetching, matrix construction) run once; only the NMF is repeated per
seed.

## Note

**Beta feature:** This function is experimental and the API may change
without notice in future versions.

## See also

[`decomposeSubnetworkByTopic`](https://vitek-lab.github.io/MSstatsBioNet/reference/decomposeSubnetworkByTopic.md)

## Examples

``` r
if (FALSE) { # \dontrun{
input <- data.table::fread(system.file(
    "extdata/groupComparisonModel.csv",
    package = "MSstatsBioNet"
))
subnetwork <- getSubnetworkFromIndra(input)
cmp <- compareTopicModels(subnetwork, seeds = 1:20, n_topics = 5)
cmp$summary
cmp$test            # p < 0.05 => PPI changes topics beyond chance
cmp$dispersion      # how stable each mode is across seeds
} # }
```
