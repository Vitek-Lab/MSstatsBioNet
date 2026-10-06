# MSstatsBioNet (development version)

## New features

* New function `validate_network()` checks a `list(nodes, edges)` network
against the edge and node contract (v1.0): required columns and types, the
statement-type and entity-type vocabularies, value ranges, `NA` statistics
on unmeasured nodes, and that every edge endpoint is a node. It stops with
an error listing every problem found.

## Breaking changes

* The edges returned by `getSubnetworkFromIndra()` follow the edge contract.
Columns are renamed, with no aliases: `stmt_hash` to `statement_id` (now
always character, taken from INDRA's full-precision `matches_hash`),
`evidenceLink` to `evidence_url`, `evidenceCount` to `evidence_count` (now
integer), and `sourceCounts` to `evidence_sources`. New columns: `directed`
(`FALSE` for symmetric types such as `Complex`), `confidence` (INDRA belief
score), `backend_database` (`"INDRA"`), and `query_type` (`"subnetwork"`).
`filterSubnetworkByContext()`, the topic-model functions, and
`cytoscapeNetwork()` read the new names, so edges built by hand for them
must use `statement_id` and `evidence_url`.
* `evidence_url` links to the INDRA page for that one statement. The old
`evidenceLink` listed every statement between the two agents.
* The `cytoscapeNetwork()` widget's `_edge_clicked` Shiny input reports
`evidence_url` in place of `evidenceLink`.
* The nodes returned by `getSubnetworkFromIndra()` follow the node contract.
Columns are renamed, with no aliases: `entityName` to `entity_name`,
`entityId` to `entity_id`, and `Site` to `site`. `logFC` keeps its name,
whatever the log base of the input. New column: `namespace`, the
grounding namespace(s) aligned with `entity_id`. The ID columns are always
character, even when every grounded ID is numeric.
* `getSubnetworkFromIndra()` calls `validate_network()` on its result, so it
stops with an error instead of returning a network that breaks the contract.
* `cytoscapeNetwork()`, `exportNetworkToHTML()`, and
`previewNetworkInBrowser()` read the new node column names (`entity_name`,
`site`). Nodes built by hand for them must use these names.
* `cytoscapeNetwork()` now errors when `displayLabelType` is not `"id"` or
`"entity_name"`. Other values used to fall back to `"id"` silently.

## Deprecated

* `cytoscapeNetwork(displayLabelType = "entityName")` (and the same argument
of `exportNetworkToHTML()` and `previewNetworkInBrowser()`) is deprecated. It
gives a warning and is treated as `"entity_name"`, the renamed nodes column.
It will become defunct in the next Bioconductor release.
* In `getSubnetworkFromIndra()`, the arguments `paper_count_cutoff`,
`correlation_cutoff`, and `protein_level_data` are deprecated. Supplying any
of them gives a warning. They will become defunct in the next Bioconductor
release and be removed in the one after.
    * `paper_count_cutoff` is now ignored. INDRA does not return paper counts
    (the `paperCount` column is always 1), so the default had no effect and
    values of 2 or more removed every edge and raised an error.
    * `correlation_cutoff` and `protein_level_data` still work during the
    deprecation period. The `paperCount` and `correlation` columns of the
    returned edges will be removed when the arguments become defunct.

## Internal changes

* `getSubnetworkFromIndra()` now runs through an internal INDRA backend
object and the S4 generic `get_network()`, the first step toward supporting
network databases other than INDRA. Its output is unchanged. The new
functions are not exported yet.
* The error for a non-character `sources_filter` in
`getSubnetworkFromIndra()` now reads "evidence_sources must be a character
vector", the name of the argument in the new API.
* Added internal functions for the entity table that the new API takes as
input: `prepare_entities()` (one row per analyte, with its entity type,
identifier system, and organism; copies `log2FC`, `log10FC`, or `logFC` to
`logFC`; stops when the input has several comparisons in `Label` and
`label` doesn't name one), `parse_ptm_sites()`, `build_grounding_table()`, and
`select_entities()` (flags rows that pass the cutoffs and drops none).
Nothing calls them yet.
* `annotateProteinInfoFromIndra()` now runs through two internal generics
of the new API: `convert_ids()`, which grounds an entity table through the
INDRA backend (CoGEx for UniProt IDs and mnemonics, Gilda for gene symbols
and chemical names), and `get_entity_properties()`, which adds the
`is_transcription_factor`, `is_kinase`, and `is_phosphatase` columns. Its
output is unchanged. The INDRA backend now also holds the Gilda URL, and
the organism of the entity table is passed to Gilda in place of a
hard-coded human taxon ID.

# MSstatsBioNet 0.99.0

* Added function `getSubnetworkFromIndra` to extract biomolecular subnetworks 
from the INDRA database.
* Added function `visualizeNetworks` to visualize networks in Cytoscape desktop.
* Added vignette

# MSstatsBioNet 0.1.0

* Added a `NEWS.md` file to track changes to the package.