# Filter a subnetwork by contextual relevance

Fetches PubMed abstracts for evidence PMIDs, scores each abstract
against a user-supplied query, and returns only the nodes, edges, and
evidence rows whose abstracts meet the scoring cutoff.

## Usage

``` r
filterSubnetworkByContext(
  nodes,
  edges,
  query,
  cutoff = NULL,
  method = c("tag_count", "cosine")
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
  `method = "cosine"`: a single character string.

- cutoff:

  Numeric threshold applied to the chosen scoring method.

  - `"tag_count"`: integer \>= 0; abstracts must contain at least this
    many tags. Max possible value is `length(query)`. Default `1`.

  - `"cosine"`: numeric in `[-1, 1]`; abstracts must score \>= this
    value. Default `0.10`.

- method:

  One of `"tag_count"` (default) or `"cosine"`.

## Value

A named list with three elements:

- nodes:

  Filtered nodes dataframe (only nodes present in kept edges)

- edges:

  Filtered edges dataframe

- evidence:

  Dataframe with columns: source, target, interaction, site,
  evidenceLink, stmt_hash, text, pmid, score. The `score` column
  contains tag counts (integer) or cosine similarities (numeric)
  depending on the method used.

## Details

Two scoring methods are available, controlled by the `method` argument:

- `"tag_count"` (default):

  Counts how many tags from `query` appear as substrings in the abstract
  (case-insensitive). The score for each abstract is an integer in
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
