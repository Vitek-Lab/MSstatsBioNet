# Filter a subnetwork by contextual relevance

Fetches PubMed abstracts for evidence PMIDs, scores each abstract
against a user-supplied query, and returns only the nodes, edges, and
evidence rows whose abstracts meet the scoring cutoff.

## Usage

``` r
filterSubnetworkByContext(
  nodes,
  edges,
  query = NULL,
  cutoff = NULL,
  method = c("tag_count", "cosine"),
  exclude_keywords = NULL
)
```

## Arguments

- nodes:

  A dataframe of network nodes.

- edges:

  A dataframe of network edges with columns: source, target,
  interaction, site, evidenceLink, stmt_hash.

- query:

  For `method = "tag_count"`: a character vector of tags, e.g.
  `c("CHEK1", "DNA damage", "DNA damage repair")`. For
  `method = "cosine"`: a single character string. May be `NULL`
  (default) when `exclude_keywords` is supplied; abstracts are then not
  scored (`score` is `NA`) and only the keyword exclusion is applied.

- cutoff:

  Numeric threshold applied to the chosen scoring method.

  - `"tag_count"`: integer \>= 0; abstracts must contain at least this
    many tags. Max possible value is `length(query)`. Default `1`.

  - `"cosine"`: numeric in `[-1, 1]`; abstracts must score \>= this
    value. Default `0.10`.

- method:

  One of `"tag_count"` (default) or `"cosine"`.

- exclude_keywords:

  Optional character vector of keywords. Abstracts containing any of
  them as a whole word or phrase (case-insensitive) are removed,
  regardless of their score. For example, `"colon"` matches "colon" and
  "colon-specific" but not "colonize" or "colons"; list variants such as
  plurals explicitly. To exclude by keyword only, omit `query`. Default
  `NULL` excludes nothing.

## Value

A named list with four elements:

- nodes:

  Filtered nodes dataframe (only nodes present in kept edges)

- edges:

  Filtered edges dataframe

- evidence:

  Dataframe with columns: source, target, interaction, site,
  evidenceLink, stmt_hash, text, pmid, score. The `score` column
  contains tag counts (integer) or cosine similarities (numeric)
  depending on the method used.

- abstracts:

  Named character vector mapping each PMID in `evidence` to its abstract
  text.

The `evidence` and `abstracts` elements can be passed to the same-named
arguments of
[`decomposeSubnetworkByTopic`](https://vitek-lab.github.io/MSstatsBioNet/reference/decomposeSubnetworkByTopic.md)
or
[`decomposeSubnetworkIntoHierarchicalTopics`](https://vitek-lab.github.io/MSstatsBioNet/reference/decomposeSubnetworkIntoHierarchicalTopics.md),
together with the returned list as `subnetwork`, so INDRA and PubMed are
not queried again.

## Details

Two scoring methods are available, controlled by the `method` argument:

- `"tag_count"` (default):

  Counts how many tags from `query` appear as whole words or phrases in
  the abstract (case-insensitive), so `"colon"` does not match "colony"
  or "colonize". The score for each abstract is an integer in
  `[0, length(query)]`. Set `cutoff` to the minimum number of tags that
  must appear - e.g. `cutoff = 2` keeps abstracts that mention at least
  2 of your tags. `query` must be a character *vector* of tags when
  using this method.

- `"cosine"`:

  Scores abstracts using TF-IDF cosine similarity against `query`.
  Scores are in `[-1, 1]` (in practice `[0, 1]` for text). Set `cutoff`
  to a decimal threshold - e.g. `cutoff = 0.10`. `query` should be a
  single character string; expand it with synonyms and related terms for
  better recall under exact token matching.

## Note

**Beta feature:** This function is experimental and the API may change
without notice in future versions.

## Examples

``` r
if (FALSE) { # \dontrun{
filtered <- filterSubnetworkByContext(
    subnetwork$nodes, subnetwork$edges,
    query = c("DNA damage", "DNA repair"),
    exclude_keywords = c("review")
)
hierarchy <- decomposeSubnetworkIntoHierarchicalTopics(
    filtered,
    evidence  = filtered$evidence,
    abstracts = filtered$abstracts
)
} # }
```
