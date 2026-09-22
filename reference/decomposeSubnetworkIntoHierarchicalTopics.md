# Recursively decompose a subnetwork into a hierarchy of topic subnetworks

Repeatedly applies
[`decomposeSubnetworkByTopic`](https://vitek-lab.github.io/MSstatsBioNet/reference/decomposeSubnetworkByTopic.md)
to its own topic subnetworks until every branch is small enough to
inspect by hand (at most `max_edges` edges), producing a topic tree:
broad themes near the root and increasingly specific sub-themes towards
the leaves.

## Usage

``` r
decomposeSubnetworkIntoHierarchicalTopics(
  subnetwork,
  max_edges = 10,
  n_topics = 5,
  edge_topic_cutoff = 0.9,
  max_depth = 5,
  evidence = NULL,
  abstracts = NULL,
  ...
)
```

## Arguments

- subnetwork:

  list with `nodes` and `edges` data.frames, e.g. the output of
  [`getSubnetworkFromIndra`](https://vitek-lab.github.io/MSstatsBioNet/reference/getSubnetworkFromIndra.md).

- max_edges:

  a branch with at most this many edges is not split further. Default
  10.

- n_topics:

  number of topics per split. Default 5.

- edge_topic_cutoff:

  topic-share threshold for assigning an edge to a topic at each split;
  see
  [`decomposeSubnetworkByTopic`](https://vitek-lab.github.io/MSstatsBioNet/reference/decomposeSubnetworkByTopic.md).
  The default of 0.9 gives a near-partition, so sibling topics rarely
  share edges. Lower values allow an edge to appear in several sibling
  topics.

- max_depth:

  maximum depth of the tree (the root is depth 0). Default 5.

- evidence:

  optional pre-fetched evidence data.frame, e.g.
  `attr(topics, "corpus")$evidence` from
  [`decomposeSubnetworkByTopic`](https://vitek-lab.github.io/MSstatsBioNet/reference/decomposeSubnetworkByTopic.md)
  or `result$corpus$evidence` from a previous call of this function.
  Default `NULL` queries INDRA once.

- abstracts:

  optional named character vector mapping PMID to abstract text. Only
  missing PMIDs are fetched from PubMed. Default `NULL`.

- ...:

  further arguments passed to
  [`decomposeSubnetworkByTopic`](https://vitek-lab.github.io/MSstatsBioNet/reference/decomposeSubnetworkByTopic.md),
  e.g. `n_top_terms`, `min_term_count`, `include_ppi`, `seed`.

## Value

An object of class `topicHierarchy`: a list with

- tree:

  data.frame with one row per topic, in depth-first order. Columns: `id`
  (`"root"`, `"1"`, `"1.2"`, ...), `parent_id`, `depth`, `topic` (index
  within the parent's split), `n_edges`, `n_nodes`, `n_papers` (papers
  supporting the topic's edges), `n_children`, `is_leaf`, `stop_reason`,
  `mean_topic_weight` (mean share of the topic's edges' loading, a
  cohesion score), `top_terms` (collapsed with `", "`), `label`, and
  `pathString` (`"root/1/1.2"`, for `data.tree`).

- subnetworks:

  named list keyed by `id`; each element is a subnetwork (`nodes`,
  `edges`, `topTerms`, `pmids`) that can be passed to
  [`cytoscapeNetwork`](https://vitek-lab.github.io/MSstatsBioNet/reference/cytoscapeNetwork.md)
  or
  [`exportNetworkToHTML`](https://vitek-lab.github.io/MSstatsBioNet/reference/exportNetworkToHTML.md).

- edge_membership:

  long data.frame with one row per (topic, edge): `id`, `depth`,
  `is_leaf`, `source`, `target`, `interaction`, `topicWeight`. Filter on
  `is_leaf` to see which fine-grained topic(s) each edge ends up in.

- corpus:

  list with the `evidence` and `abstracts` used, for reuse via the
  `evidence` and `abstracts` arguments.

- params:

  the settings used.

## Details

INDRA evidence and PubMed abstracts are gathered once for the input
subnetwork and reused for every sub-decomposition, so no further network
requests are made during the recursion. Each sub-decomposition rebuilds
its vocabulary and refits the NMF on only the papers supporting that
branch's edges, which lets finer topics emerge.

A branch stops splitting (becomes a leaf) when any of the following
holds, recorded in the `stop_reason` column of `tree`:

- small_enough:

  it has at most `max_edges` edges.

- max_depth:

  it sits at depth `max_depth`.

- too_few_papers:

  fewer than two papers support its edges.

- no_split:

  every child topic contained all of its edges (or none), so splitting
  would make no progress.

- failed: \<message\>:

  the decomposition raised an error, e.g. no usable words in the
  branch's abstracts.

Edges without any PMID-backed evidence cannot be assigned to a topic and
only appear at the root.

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
hierarchy <- decomposeSubnetworkIntoHierarchicalTopics(
    subnetwork, max_edges = 10, n_topics = 5, edge_topic_cutoff = 0.9
)
hierarchy                       # indented topic tree
leaves <- hierarchy$tree[hierarchy$tree$is_leaf, ]
leaves[order(-leaves$mean_topic_weight), c("id", "n_edges", "top_terms")]

# Inspect one fine-grained topic as a network.
leaf <- hierarchy$subnetworks[["1.2"]]
exportNetworkToHTML(leaf$nodes, leaf$edges)

# Tree visualization with other packages, e.g.
# data.tree::as.Node(hierarchy$tree)
# igraph::graph_from_data_frame(
#     hierarchy$tree[-1, c("parent_id", "id")])
} # }
```
