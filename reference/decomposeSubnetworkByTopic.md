# Decompose a subnetwork into topic-specific subnetworks via joint NMF

Takes a subnetwork (the output of
[`getSubnetworkFromIndra`](https://vitek-lab.github.io/MSstatsBioNet/reference/getSubnetworkFromIndra.md))
and splits it into a list of smaller, topic-specific subnetworks
discovered with unsupervised non-negative matrix factorization (NMF).

## Usage

``` r
decomposeSubnetworkByTopic(
  subnetwork,
  n_topics = 5,
  edge_topic_cutoff = 0.2,
  n_top_terms = 10,
  min_term_count = 2,
  max_iter = 200,
  tol = 1e-04,
  seed = 1,
  include_ppi = TRUE
)
```

## Arguments

- subnetwork:

  list with `nodes` and `edges` data.frames, e.g. the output of
  [`getSubnetworkFromIndra`](https://vitek-lab.github.io/MSstatsBioNet/reference/getSubnetworkFromIndra.md).

- n_topics:

  number of topics (rank of the factorization). Default 5.

- edge_topic_cutoff:

  numeric in `[0, 1]`; an edge is added to a topic's subnetwork when the
  topic carries at least this share of the edge's total loading. Each
  edge is always included in at least its highest-loading topic. Default
  0.2.

- n_top_terms:

  number of top words to report per topic. Default 10.

- min_term_count:

  minimum corpus frequency for a word to be kept when building `X_text`.
  Default 2.

- max_iter:

  maximum number of NMF multiplicative-update iterations. Default 200.

- tol:

  relative-change tolerance for NMF early stopping. Default 1e-4.

- seed:

  random seed for NMF initialization. Default 1.

- include_ppi:

  logical; if `TRUE` (default) the PPI/edge matrix is factorized jointly
  with the text matrix via a shared basis. If `FALSE`, NMF is run on the
  paper-word matrix only and edge-topic loadings are derived afterwards
  by folding edge counts onto the text-learned topics, so the PPIs do
  not influence the topics themselves.

## Value

A list of length `n_topics`, named `topic_1` ... `topic_k`. Each element
is a topic-specific subnetwork: a list with

- nodes:

  nodes data.frame restricted to the topic's edges.

- edges:

  edges data.frame for the topic, with an added `topicWeight` column
  (the edge's topic share).

- topic:

  the topic index.

- topTerms:

  character vector of the topic's top words.

- pmids:

  PMIDs whose strongest topic loading is this topic.

The full factorization (W, H_text, H_edges, etc.) is attached as the
`"nmf"` attribute of the returned list.

## Details

The procedure is:

1.  For every edge, the supporting INDRA evidence is retrieved and the
    PubMed abstract of each referenced PMID is fetched. Papers (PMIDs)
    are the shared unit of analysis.

2.  Two matrices are built that share the same rows (papers): `X_text`
    (papers x words, term counts from the abstracts) and `X_edges`
    (papers x unique `source_target_interaction` combinations,
    evidence-sentence counts).

3.  NMF learns a basis matrix `W` (papers x topics). When
    `include_ppi = TRUE` (the default) a *joint* NMF learns a single
    shared `W` such that \\X\_{text} \approx W H\_{text}\\ and
    \\X\_{edges} \approx W H\_{edges}\\, tying each learned topic to
    both a set of words and a set of edges. When `include_ppi = FALSE`
    the factorization uses only `X_text` (\\X\_{text} \approx W
    H\_{text}\\); the PPI evidence is excluded from the modeling and
    edge-topic loadings are instead derived afterwards by folding the
    edge counts onto the text-learned topics (\\H\_{edges} = W^\top
    X\_{edges}\\). This lets you compare topic structure with and
    without the PPI view.

4.  Each topic becomes its own subnetwork: an edge is included in a
    topic when that topic carries at least `edge_topic_cutoff` of the
    edge's loading (soft, overlapping assignment), and nodes are
    restricted to those touched by the kept edges.

## Note

**Beta feature:** This function is experimental and the API may change
without notice in future versions.

## See also

[`getSubnetworkFromIndra`](https://vitek-lab.github.io/MSstatsBioNet/reference/getSubnetworkFromIndra.md),
[`filterSubnetworkByContext`](https://vitek-lab.github.io/MSstatsBioNet/reference/filterSubnetworkByContext.md)

## Examples

``` r
if (FALSE) { # \dontrun{
input <- data.table::fread(system.file(
    "extdata/groupComparisonModel.csv",
    package = "MSstatsBioNet"
))
subnetwork <- getSubnetworkFromIndra(input)
topics <- decomposeSubnetworkByTopic(subnetwork, n_topics = 5)
topics$topic_1$topTerms
exportNetworkToHTML(topics$topic_1$nodes, topics$topic_1$edges)
} # }
```
