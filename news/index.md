# Changelog

## MSstatsBioNet (development version)

### New features

- New function
  [`validate_network()`](https://vitek-lab.github.io/MSstatsBioNet/reference/validate_network.md)
  checks a `list(nodes, edges)` network against the edge and node
  contract (v1.0): required columns and types, the statement-type and
  entity-type vocabularies, value ranges, `NA` statistics on unmeasured
  nodes, and that every edge endpoint is a node. It stops with an error
  listing every problem found.

### Deprecated

- In
  [`getSubnetworkFromIndra()`](https://vitek-lab.github.io/MSstatsBioNet/reference/getSubnetworkFromIndra.md),
  the arguments `paper_count_cutoff`, `correlation_cutoff`, and
  `protein_level_data` are deprecated. Supplying any of them gives a
  warning. They will become defunct in the next Bioconductor release and
  be removed in the one after.
  - `paper_count_cutoff` is now ignored. INDRA does not return paper
    counts (the `paperCount` column is always 1), so the default had no
    effect and values of 2 or more removed every edge and raised an
    error.
  - `correlation_cutoff` and `protein_level_data` still work during the
    deprecation period. The `paperCount` and `correlation` columns of
    the returned edges will be removed when the arguments become
    defunct.

## MSstatsBioNet 0.99.0

- Added function `getSubnetworkFromIndra` to extract biomolecular
  subnetworks from the INDRA database.
- Added function `visualizeNetworks` to visualize networks in Cytoscape
  desktop.
- Added vignette

## MSstatsBioNet 0.1.0

- Added a `NEWS.md` file to track changes to the package.
