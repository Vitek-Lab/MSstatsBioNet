# MSstatsBioNet (development version)

## New features

* A new API for building networks from MSstats results, which separates
the steps that `annotateProteinInfoFromIndra()` and
`getSubnetworkFromIndra()` combine, and works for proteins, PTM sites, and
metabolites:
    * `prepare_entities()` builds an entity table with one row per analyte,
    its entity type (`"protein"`, `"ptm_site"`, `"metabolite"`, ...),
    identifier system, organism, and statistics. It copies `log2FC`,
    `log10FC`, or `logFC` to `logFC`, parses PTM sites, and stops when the
    input has several comparisons in `Label` unless `label` names one.
    * `indra_backend()` creates the INDRA backend, and
    `backend_capabilities()` lists what a backend supports.
    * `convert_ids()` grounds the entity table (CoGEx for UniProt IDs and
    mnemonics, Gilda for gene symbols and chemical names), and
    `get_entity_properties()` adds `is_transcription_factor`, `is_kinase`,
    and `is_phosphatase`.
    * `select_entities()` flags the rows that pass the cutoffs, and drops
    none, so that nodes in the input are recognized even when they fail
    the cutoffs.
    * `get_network(backend, entities, query = subnetwork_query())` queries
    the backend and returns `nodes` and `edges` that meet the contract. It
    takes `interaction_types`, `min_evidence`, `min_confidence` (new: drops
    edges below it, and edges with no score), `evidence_sources`, and
    `include_entities`. It prints the question it asks as a message, e.g.
    "INDRA subnetwork: how are 42 selected proteins connected to each
    other, with no other nodes added?".
    * `?network_queries` describes the questions a query can ask, with a
    glossary. `subnetwork_query()` is the first; more are planned.
    * The S4 classes `NetworkBackend`, `IndraBackend`, `NetworkQuery`, and
    `SubnetworkQuery` are exported, so other packages can add backends.
* New function `validate_network()` checks a `list(nodes, edges)` network
against the edge and node contract (v1.0): required columns and types, the
statement-type and entity-type vocabularies, value ranges, `NA` statistics
on nodes that are not in the input data, and that every edge's source and
target are nodes. It stops with an error listing every problem found.

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
* The nodes returned by `getSubnetworkFromIndra()` say whether each node
is in the input data and why it is in the network. New columns: `entity_type`
(`"protein"`, `"ptm_site"`, `"metabolite"`, `"family"`, ...), `measured`
(`TRUE` for analytes in `input`), `included_in_query`, `node_role`
(`"passed_cutoffs"`, or `"user_added"` for nodes there only through
`force_include_other`), and `has_measured_sites` (`TRUE` for proteins with
PTM site rows). `nodes` is now always a data.frame, also when `input` is a
data.table.
* Nodes added through `force_include_other` that are not in `input` have
`NA` for `logFC` and `adj.pvalue`, and `measured = FALSE`. They used to get
`logFC = 0` and `adj.pvalue = 1`, which looked like a protein in the input data
with no change. They now also carry INDRA's namespace and identifier, and an
`entity_type` from the namespace (FamPlex families are `"family"`).
* INDRA statements between identifiers that are equal but in different
namespaces (`HGNC:1234` and `CHEBI:1234`) are no longer merged into one
edge. The order of the edge rows can differ from earlier versions.
* When an INDRA node matches rows of `input` that belong to different
nodes, a message names it. It keeps INDRA's name as its `id`, as before,
and now has `NA` statistics.
* `getSubnetworkFromIndra()` stops when `input` has more than one
comparison in its `Label` column, or repeats a `Protein` value. Filter
`input` to one comparison first.
* PTM sites are parsed from the end of each `;`-separated member of
`Protein`: `CLH1_HUMAN_S148` gives site `S148` (was `HUMAN_S148`),
`P1_S148;P2_T5` gives sites `S148;T5` (was `S148;P2_T5`), and an
identifier whose last part is not a site (`P1_S148_extra`) has none. A PTM
row without a `GlobalProtein` column is drawn on its parent protein's node
(`P1`), not on a node of its own (`P1_S148`).
* The message about rows with no grounding now counts only the rows that
pass the cutoffs, since only those are sent to INDRA.
* `validate_network()` requires the node columns `entity_type`,
`entity_name`, `namespace`, `entity_id`, `measured`, `included_in_query`,
and `node_role`.

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
network databases other than INDRA.
* The error for a non-character `sources_filter` in
`getSubnetworkFromIndra()` now reads "evidence_sources must be a character
vector", the name of the argument in the new API.
* `annotateProteinInfoFromIndra()` now runs through `convert_ids()` and
`get_entity_properties()`. Its output is unchanged. The INDRA backend now also holds the Gilda URL, and
the organism of the entity table is passed to Gilda in place of a
hard-coded human taxon ID.
* `getSubnetworkFromIndra()` builds an entity table from `input`, flags
the rows to query with `select_entities()`, and passes all rows to
`get_network()`, which matches the nodes INDRA returns against every row.
Like `get_network()`, it now prints the question it asks as a message.

# MSstatsBioNet 0.99.0

* Added function `getSubnetworkFromIndra` to extract biomolecular subnetworks 
from the INDRA database.
* Added function `visualizeNetworks` to visualize networks in Cytoscape desktop.
* Added vignette

# MSstatsBioNet 0.1.0

* Added a `NEWS.md` file to track changes to the package.