# Filters that run on a finished network, after get_network().

#' Remove evidence curated as incorrect
#'
#' Subtracts the evidence curated as incorrect, from
#' \code{\link{get_curations}()}, from each edge's \code{evidence_count},
#' then drops the edges left below \code{min_evidence}, and the nodes those
#' edges leave without any edge. Edges left with no evidence are dropped
#' whatever \code{min_evidence} is. A message says how many edges were
#' dropped.
#'
#' Run it after any filter that pools edges, such as the PTM-site filter of
#' \code{\link{getSubnetworkFromIndra}()}: dropping edges first can change
#' which edges such a filter keeps.
#'
#' @param network list of \code{nodes} and \code{edges} from
#' \code{\link{get_network}()}. \code{edges} needs the columns
#' \code{source}, \code{target}, \code{statement_id}, and
#' \code{evidence_count}, and \code{backend_database} unless
#' \code{backend} is given.
#' @param min_evidence minimum evidence count per edge, after subtracting
#' the incorrect evidence
#' @param backend the backend to get the curations from, e.g.
#' \code{\link{indra_backend}()}. \code{NULL} (default) uses the default
#' backend named in each edge's \code{backend_database}. Pass a backend to
#' use one built with non-default settings, such as
#' \code{indra_backend(curation_url = )}. A backend without curations, such
#' as one with no \code{get_curations()} method, gives an error naming it.
#' @return \code{network}, with \code{edges$evidence_count} reduced and the
#' dropped edges and nodes removed. Other elements of \code{network} are
#' kept.
#' @seealso \code{\link{get_curations}()}
#' @export
#' @examples
#' \donttest{
#' input <- data.table::fread(system.file(
#'     "extdata/groupComparisonModel.csv",
#'     package = "MSstatsBioNet"
#' ))
#' indra <- indra_backend()
#' entities <- prepare_entities(input, entity_type = "protein",
#'                              id_type = "uniprot")
#' entities <- convert_ids(indra, entities)
#' entities <- select_entities(entities, pvalue_cutoff = 0.05)
#' network <- get_network(indra, entities, interaction_types = "Complex")
#' network <- filter_by_curation(network, min_evidence = 2)
#' head(network$edges)
#' }
filter_by_curation <- function(network, min_evidence = 1, backend = NULL) {
    .check_filter_by_curation_input(network, min_evidence)
    .check_backend_argument(backend)
    edges <- network$edges
    incorrect_counts <- .count_incorrect_evidence(edges, backend)
    edges$evidence_count <- as.integer(edges$evidence_count - incorrect_counts)
    keep <- edges$evidence_count >= min_evidence & edges$evidence_count >= 1
    if (any(!keep)) {
        message("Dropping ", sum(!keep), " edge(s) with fewer than ",
                max(min_evidence, 1), " evidence after removing the ",
                "evidence curated as incorrect.")
    }
    .keep_edges(network, keep, edges)
}

#' Keep some edges of a network, and drop the nodes they leave without edges
#'
#' Nodes that had no edge before are kept, as are the other elements of
#' the network.
#' @param network network list
#' @param keep logical, one per edge
#' @param edges the edges to subset, \code{network$edges} by default
#' @return \code{network}
#' @keywords internal
#' @noRd
.keep_edges <- function(network, keep, edges = network$edges) {
    nodes <- network$nodes
    endpoints_before <- c(edges$source, edges$target)
    edges <- edges[keep, , drop = FALSE]
    endpoints_after <- c(edges$source, edges$target)
    left_without_edges <- nodes$id %in% endpoints_before &
        !nodes$id %in% endpoints_after
    network$nodes <- nodes[!left_without_edges, , drop = FALSE]
    network$edges <- edges
    network
}

#' Check the input of filter_by_curation()
#' @param network the network argument
#' @param min_evidence the min_evidence argument
#' @return \code{NULL}, invisibly; errors otherwise
#' @keywords internal
#' @noRd
.check_filter_by_curation_input <- function(network, min_evidence) {
    if (!is.list(network) || !is.data.frame(network$nodes) ||
        !is.data.frame(network$edges)) {
        stop("`network` must be a list with `nodes` and `edges` ",
             "data.frames, e.g. from get_network().", call. = FALSE)
    }
    missing_cols <- setdiff(c("source", "target", "statement_id",
                              "evidence_count"), names(network$edges))
    if (length(missing_cols) > 0) {
        stop("`network$edges` is missing the column(s): ",
             paste(missing_cols, collapse = ", "), call. = FALSE)
    }
    if (!"id" %in% names(network$nodes)) {
        stop("`network$nodes` is missing the column: id", call. = FALSE)
    }
    if (!is.numeric(min_evidence) || length(min_evidence) != 1 ||
        !is.finite(min_evidence)) {
        stop("`min_evidence` must be a single number.", call. = FALSE)
    }
    invisible(NULL)
}

#' Keep the edges that modify a PTM site in the network
#'
#' Keeps an edge when its \code{site} is one of the PTM sites of its
#' target, that is, of a node row with the target's \code{id}. For a
#' network built from MSstatsPTM results, these are the edges that explain
#' a site in the input, e.g. a kinase phosphorylating a site whose
#' phosphorylation changed. Each edge is checked on its own, so the order
#' of this filter and \code{\link{filter_by_curation}()} doesn't matter.
#'
#' Edges without a \code{site} are dropped. These include every edge that
#' isn't a modification (\code{Complex}, \code{Activation}, ...), and every
#' edge from a backend that reports no sites, such as STRING: filter an
#' INDRA network before merging it with one from such a backend, if those
#' edges should stay. Nodes are dropped when they had edges and lost them
#' all. Nodes that had no edge, and other elements of \code{network}, such
#' as \code{provenance}, are kept.
#'
#' Node sites are split on \code{"_"} (several sites of one peptide, e.g.
#' \code{"S1039_S1042"}) and \code{";"} (members of a protein group);
#' edge sites on \code{";"}.
#'
#' @param network list of \code{nodes} and \code{edges} from
#' \code{\link{get_network}()}, with PTM site rows in \code{nodes}
#' (entity tables with \code{entity_type = "ptm_site"})
#' @return \code{network}, with the edges whose site is not a site of
#' their target removed, and the nodes they leave without edges. Stops
#' when no edges remain.
#' @seealso \code{\link{add_entity_values}()}
#' @export
#' @examples
#' \donttest{
#' input <- data.table::fread(system.file(
#'     "extdata/garrido-2024.csv",
#'     package = "MSstatsBioNet"
#' ))
#' indra <- indra_backend()
#' entities <- prepare_entities(input, entity_type = "ptm_site",
#'                              id_type = "uniprot")
#' entities <- convert_ids(indra, entities)
#' entities <- select_entities(entities, pvalue_cutoff = 0.05,
#'                             logfc_cutoff = 1)
#' network <- get_network(indra, entities,
#'                        interaction_types = "Phosphorylation",
#'                        include_entities = "HGNC:3236")
#' network <- filter_by_ptm_site(network)
#' head(network$edges)
#' }
filter_by_ptm_site <- function(network) {
    validate_network(network)
    nodes <- network$nodes
    node_sites <- .get_column_or_default(nodes, "site", NA_character_)
    site_rows <- !is.na(node_sites)
    if (!any(site_rows)) {
        stop("network$nodes has no PTM site rows (site), so ",
             "filter_by_ptm_site() would drop every edge. Build the network ",
             "from an entity table with entity_type = \"ptm_site\".",
             call. = FALSE)
    }
    sites_by_node <- lapply(
        split(strsplit(node_sites[site_rows], "[_;]"), nodes$id[site_rows]),
        function(sites) unique(trimws(unlist(sites))))
    edges <- network$edges
    edge_sites <- .get_column_or_default(edges, "site", NA_character_)
    keep <- vapply(seq_len(nrow(edges)), function(i) {
        if (is.na(edge_sites[i])) {
            return(FALSE)
        }
        sites <- trimws(strsplit(edge_sites[i], ";", fixed = TRUE)[[1]])
        any(sites[nzchar(sites)] %in% sites_by_node[[edges$target[i]]])
    }, logical(1))
    if (!any(keep)) {
        stop("No edges remain after filtering by PTM site: no edge's site ",
             "is a site of its target in the network.", call. = FALSE)
    }
    if (any(!keep)) {
        message("Dropping ", sum(!keep), " edge(s) whose site is not a PTM ",
                "site of their target in the network.")
    }
    .keep_edges(network, keep)
}

#' Copy the values of an entity table onto a network's nodes
#'
#' Sets \code{logFC} and \code{adj.pvalue} of the nodes that match a row
#' of \code{entities}, without querying a backend. Use it to show the
#' results of another experiment on a network: for example, a network built
#' from MSstatsPTM's \code{ADJUSTED.Model}, with its protein nodes colored
#' by the global profiling results (\code{PROTEIN.Model}).
#'
#' A node row matches an entity row with the same node \code{id} and
#' \code{site} (for a PTM site row, its parent protein and site; for a
#' protein, \code{site = NA}). Node rows that match no entity row by
#' \code{id} are then matched by grounding (\code{namespace} and
#' \code{entity_id}) and \code{site}. This finds nodes whose \code{id} is
#' the backend's name, such as a kinase added with
#' \code{get_network(include_entities =)} that INDRA calls \code{"SRC"},
#' while the entity table has \code{"P12931"}. An entity row that matched a
#' node by \code{id} is not used again by grounding, so an isoform node
#' doesn't take the values of its canonical protein. A node row that
#' matches several entity rows gets no values, with a message.
#'
#' A protein that has only PTM site rows in the network, and matches a
#' protein-level entity row, gets a new protein-level row (\code{site =
#' NA}) with that row's values, so \code{\link{cytoscapeNetwork}()} colors
#' the protein node by it.
#'
#' Matched node rows get \code{measured = TRUE}. Rows of \code{entities}
#' that match no node are not added: \code{\link{get_network}()} decides
#' which nodes are in the network. Edges and the other elements of
#' \code{network} are unchanged. The node table does not record which
#' experiment the values come from.
#'
#' @param network list of \code{nodes} and \code{edges} from
#' \code{\link{get_network}()}
#' @param entities entity table from \code{\link{prepare_entities}()},
#' grounded by \code{\link{convert_ids}()}. Its \code{included_in_query}
#' column is ignored.
#' @param overwrite whether to replace the \code{logFC} and
#' \code{adj.pvalue} of node rows that already have them. \code{FALSE}
#' (default) stops and names those nodes.
#' @return \code{network}, with the values of \code{entities} on the
#' matching nodes
#' @seealso \code{\link{filter_by_ptm_site}()}, \code{\link{merge_networks}()}
#' @export
#' @examples
#' # A PTM network: phosphorylation of two EGFR sites by SRC, which is not
#' # in the PTM results, so it has INDRA's name as its id
#' network <- list(
#'     nodes = data.frame(
#'         id = c("P00533", "P00533", "SRC"),
#'         entity_type = c("ptm_site", "ptm_site", "protein"),
#'         entity_name = c("EGFR", "EGFR", "SRC"),
#'         namespace = "HGNC", entity_id = c("3236", "3236", "11283"),
#'         measured = c(TRUE, TRUE, FALSE), included_in_query = TRUE,
#'         node_role = c("passed_cutoffs", "passed_cutoffs", "user_added"),
#'         site = c("S1039", "Y1092", NA),
#'         has_measured_sites = c(TRUE, TRUE, FALSE),
#'         logFC = c(1.5, -0.4, NA), adj.pvalue = c(0.01, 0.3, NA)),
#'     edges = data.frame(
#'         source = "SRC", target = "P00533", interaction = "Phosphorylation",
#'         directed = TRUE, site = "Y1092", confidence = 0.9,
#'         evidence_count = 3L, evidence_url = "https://example.org/1",
#'         statement_id = "1", backend_database = "INDRA",
#'         query_type = "subnetwork")
#' )
#' # Global profiling results, grounded (convert_ids() would fill these)
#' global <- prepare_entities(
#'     data.frame(Protein = c("P00533", "P12931"), log2FC = c(0.2, 1.1),
#'                adj.pvalue = c(0.6, 0.02)),
#'     entity_type = "protein", id_type = "uniprot")
#' global$namespace <- "HGNC"
#' global$entity_id <- c("3236", "11283")
#' global$entity_name <- c("EGFR", "SRC")
#' network <- add_entity_values(network, global)
#' network$nodes[, c("id", "site", "measured", "logFC")]
add_entity_values <- function(network, entities, overwrite = FALSE) {
    .check_add_entity_values_input(network, entities, overwrite)
    nodes <- as.data.frame(network$nodes, stringsAsFactors = FALSE)
    nodes$site <- .get_column_or_default(nodes, "site", NA_character_)
    protein_rows <- .build_missing_protein_rows(nodes)
    candidates <- rbind(nodes, protein_rows)
    matched <- .match_nodes_to_entities(candidates, entities)
    is_new <- seq_len(nrow(candidates)) > nrow(nodes)
    matched_existing <- !is.na(matched) & !is_new
    has_values <- !is.na(.get_column_or_default(candidates, "logFC",
                                                NA_real_)) |
        !is.na(.get_column_or_default(candidates, "adj.pvalue", NA_real_))
    if (!overwrite && any(matched_existing & has_values)) {
        stop("Node(s) ", .list_values_for_message(unique(
            candidates$id[matched_existing & has_values])),
            " already have logFC or adj.pvalue values. Use overwrite = TRUE ",
            "to replace them.", call. = FALSE)
    }
    keep <- !is_new | !is.na(matched)
    rows <- matched[keep]
    nodes <- candidates[keep, , drop = FALSE]
    has_match <- !is.na(rows)
    for (column in c("logFC", "adj.pvalue")) {
        values <- .get_column_or_default(nodes, column, NA_real_)
        values[has_match] <- as.numeric(.get_column_or_default(
            entities, column, NA_real_)[rows[has_match]])
        nodes[[column]] <- values
    }
    nodes$measured[has_match] <- TRUE
    added <- is_new[keep]
    nodes$entity_type[added] <- entities$entity_type[rows[added]]
    # Each new protein-level row goes just before the first row of its node
    position <- seq_len(nrow(nodes))
    position[added] <- match(nodes$id[added], nodes$id) - 0.5
    nodes <- nodes[order(position), , drop = FALSE]
    rownames(nodes) <- NULL
    message("Copied the values of entities to ", sum(has_match & !added),
            " node row(s), and added a protein-level row to ", sum(added),
            " protein(s) with only PTM site rows.")
    network$nodes <- nodes
    validate_network(network)
    network
}

#' Check the input of add_entity_values()
#' @param network,entities,overwrite the arguments
#' @return \code{NULL}, invisibly; errors otherwise
#' @keywords internal
#' @noRd
.check_add_entity_values_input <- function(network, entities, overwrite) {
    validate_network(network)
    .validate_entities(entities)
    if (!is.logical(overwrite) || length(overwrite) != 1 || is.na(overwrite)) {
        stop("`overwrite` must be TRUE or FALSE.", call. = FALSE)
    }
    grounded <- !is.na(entities$namespace) & !is.na(entities$entity_id)
    if (!any(grounded)) {
        stop("No row of entities has a grounding. Run convert_ids() on ",
             "entities first.", call. = FALSE)
    }
    if (any(!grounded)) {
        message(sum(!grounded), " row(s) of entities have no grounding and ",
                "can match nodes only by id.")
    }
    invisible(NULL)
}

#' Build a protein-level row for each protein with only PTM site rows
#'
#' A copy of the protein's first row with \code{site = NA}, the roles of
#' all its rows, \code{measured = FALSE}, and no statistics, until
#' add_entity_values() fills it from a matching entity row.
#' @param nodes nodes data.frame with a \code{site} column
#' @return nodes data.frame, one row per such protein
#' @keywords internal
#' @noRd
.build_missing_protein_rows <- function(nodes) {
    has_protein_row <- unique(nodes$id[is.na(nodes$site)])
    site_only_ids <- setdiff(unique(nodes$id), has_protein_row)
    rows <- nodes[match(site_only_ids, nodes$id), , drop = FALSE]
    rows$site <- rep(NA_character_, nrow(rows))
    rows$measured <- rep(FALSE, nrow(rows))
    for (column in intersect(c("logFC", "adj.pvalue"), colnames(rows))) {
        rows[[column]] <- rep(NA_real_, nrow(rows))
    }
    rows$node_role <- vapply(site_only_ids, function(id) {
        .join_unique_values(nodes$node_role[nodes$id == id])
    }, character(1), USE.NAMES = FALSE)
    rows$included_in_query <- vapply(site_only_ids, function(id) {
        .any_or_na(nodes$included_in_query[nodes$id == id])
    }, logical(1), USE.NAMES = FALSE)
    rows
}

#' Match node rows to entity rows, by id and site, then by grounding and
#' site
#'
#' @param nodes nodes data.frame with a \code{site} column
#' @param entities entity table
#' @return integer vector, one per node row: the matching entity row, or
#' \code{NA} when none or several match (with a message for several)
#' @keywords internal
#' @noRd
.match_nodes_to_entities <- function(nodes, entities) {
    entity_sites <- .get_column_or_default(entities, "site", NA_character_)
    entity_keys <- .build_node_keys(.get_node_ids(entities), entity_sites)
    rows_by_key <- split(seq_len(nrow(entities)), entity_keys)
    by_id <- rows_by_key[.build_node_keys(nodes$id, nodes$site)]
    matched <- vapply(by_id, function(rows) {
        if (length(rows) == 1) rows else NA_integer_
    }, integer(1), USE.NAMES = FALSE)
    ambiguous <- lengths(by_id) > 1

    groundings <- build_grounding_table(entities)
    grounding_rows <- match(groundings$id, entities$id)
    grounding_keys <- .build_node_keys(
        paste(groundings$namespace, groundings$entity_id, sep = ":"),
        entity_sites[grounding_rows])
    available <- !grounding_rows %in% matched
    rows_by_grounding <- split(grounding_rows[available],
                               grounding_keys[available])
    for (i in which(is.na(matched) & !ambiguous)) {
        if (is.na(nodes$namespace[i]) || is.na(nodes$entity_id[i])) {
            next
        }
        node_groundings <- paste(
            strsplit(nodes$namespace[i], ";", fixed = TRUE)[[1]],
            strsplit(nodes$entity_id[i], ";", fixed = TRUE)[[1]], sep = ":")
        rows <- unique(unlist(rows_by_grounding[
            .build_node_keys(node_groundings, nodes$site[i])]))
        if (length(rows) == 1) {
            matched[i] <- rows
        } else if (length(rows) > 1) {
            ambiguous[i] <- TRUE
        }
    }
    if (any(ambiguous)) {
        message(sum(ambiguous), " node row(s) match several rows of entities ",
                "and got no values: ",
                .list_values_for_message(unique(nodes$id[ambiguous])), ".")
    }
    matched
}

#' Combine networks from several queries or backends
#'
#' Combines the networks from several \code{\link{get_network}()} calls,
#' for example two queries against one backend, or one query against two
#' backends, into one network that meets the same contract
#' (\code{\link{validate_network}()}).
#'
#' Edges are the same edge when they have the same \code{source},
#' \code{target}, \code{backend_database}, and \code{statement_id}. One row
#' is kept per edge: the copy with the largest \code{evidence_count} (the
#' first one on ties), with the \code{query_type} values of all copies
#' joined by \code{";"}, e.g. \code{"subnetwork;mediated"}. A message names
#' the statements whose copies differ in \code{evidence_count} or
#' \code{confidence}. The same relation from two backends stays two rows,
#' one per \code{backend_database}, since each has its own
#' \code{confidence} and \code{evidence_url}; confidences are never
#' combined across backends. Columns that only some networks have are
#' filled with \code{NA}.
#'
#' Nodes are the same node when they have the same \code{id} and
#' \code{site} (a protein has one row per PTM site). Their
#' \code{node_role} values are joined by \code{";"}, e.g.
#' \code{"passed_cutoffs;mediator"}, \code{included_in_query} is
#' \code{TRUE} if it is in any network, and other columns take the first
#' non-\code{NA} value. \code{has_measured_sites} describes the protein, so
#' it is \code{TRUE} on every row of a protein that has it in any network:
#' merging a protein network with a PTM network marks the protein's
#' protein-level row too.
#' \code{included_in_query} describes the query, so it can differ between
#' networks built from one entity table: a node can be a latent regulator
#' in one query and part of another.
#'
#' Without \code{entities}, \code{measured}, \code{logFC}, and
#' \code{adj.pvalue} must agree for each node, because networks built from
#' one entity table always agree; otherwise \code{merge_networks()} stops
#' and names the nodes. To merge networks built from different entity
#' tables, pass the entity table as \code{entities}: \code{measured},
#' \code{logFC}, \code{adj.pvalue}, and \code{has_measured_sites} are then
#' recomputed from it (nodes not in it are \code{measured = FALSE} with
#' \code{NA} statistics).
#'
#' The \code{regulators} and \code{provenance} tables of the networks, when
#' present, are combined by row, with columns filled with \code{NA}. Other
#' list elements are dropped, with a message.
#'
#' @param ... networks: lists of \code{nodes} and \code{edges} from
#' \code{\link{get_network}()}
#' @param entities \code{NULL} (default), or the entity table from
#' \code{\link{prepare_entities}()} to recompute node status from
#' @return list of \code{nodes} and \code{edges}, plus \code{regulators}
#' and \code{provenance} when any network has them
#' @seealso \code{\link{get_network}()}, \code{\link{validate_network}()}
#' @export
#' @examples
#' \donttest{
#' input <- data.table::fread(system.file(
#'     "extdata/groupComparisonModel.csv",
#'     package = "MSstatsBioNet"
#' ))
#' indra <- indra_backend()
#' entities <- prepare_entities(input, entity_type = "protein",
#'                              id_type = "uniprot")
#' entities <- convert_ids(indra, entities)
#' # Two selections from one entity table: the edges among the proteins
#' # that pass both cutoffs are in both networks
#' loose <- get_network(indra, select_entities(entities, pvalue_cutoff = 0.05))
#' strict <- get_network(indra, select_entities(entities, pvalue_cutoff = 0.01))
#' network <- merge_networks(loose, strict)
#' # The shared edges come back once
#' c(loose = nrow(loose$edges), strict = nrow(strict$edges),
#'   merged = nrow(network$edges))
#' }
merge_networks <- function(..., entities = NULL) {
    networks <- list(...)
    .check_merge_networks_input(networks, entities)
    nodes <- .merge_nodes(lapply(networks, `[[`, "nodes"), entities)
    network <- list(nodes = nodes,
                    edges = .merge_edges(lapply(networks, `[[`, "edges")))
    for (element in MERGED_TABLE_ELEMENTS) {
        tables <- Filter(Negate(is.null), lapply(networks, `[[`, element))
        if (length(tables) > 0) {
            network[[element]] <- .bind_rows_with_fill(tables)
        }
    }
    dropped <- setdiff(unlist(lapply(networks, names)),
                       c("nodes", "edges", MERGED_TABLE_ELEMENTS))
    if (length(dropped) > 0) {
        message("merge_networks() dropped the network element(s): ",
                .list_values_for_message(unique(dropped)), ".")
    }
    validate_network(network)
    network
}

#' List elements of a network, besides nodes and edges, that
#' merge_networks() combines by row
#' @keywords internal
#' @noRd
MERGED_TABLE_ELEMENTS <- c("regulators", "provenance")

#' Node columns that must agree across networks built from one entity table
#'
#' Not included_in_query, which describes the query rather than the entity
#' table.
#' @keywords internal
#' @noRd
NODE_STATUS_COLUMNS <- c("measured", "logFC", "adj.pvalue")

#' Check the input of merge_networks()
#'
#' Each network must meet the contract on its own, so that an error names
#' the network it comes from.
#' @param networks list of the networks passed in \code{...}
#' @param entities the entities argument
#' @return \code{NULL}, invisibly; errors otherwise
#' @keywords internal
#' @noRd
.check_merge_networks_input <- function(networks, entities) {
    if (length(networks) == 0) {
        stop("merge_networks() needs at least one network.", call. = FALSE)
    }
    labels <- names(networks)
    if (is.null(labels)) {
        labels <- rep("", length(networks))
    }
    labels <- ifelse(labels == "", paste0("network ", seq_along(networks)),
                     paste0("network '", labels, "'"))
    for (i in seq_along(networks)) {
        tryCatch(validate_network(networks[[i]]), error = function(e) {
            stop("In merge_networks(), ", labels[i], ": ",
                 conditionMessage(e), call. = FALSE)
        })
    }
    if (!is.null(entities)) {
        .validate_entities(entities)
    }
    invisible(NULL)
}

#' Combine and deduplicate edges
#' @param edge_tables list of edges data.frames
#' @return edges data.frame, one row per source, target,
#' backend_database, and statement_id
#' @keywords internal
#' @noRd
.merge_edges <- function(edge_tables) {
    edges <- .bind_rows_with_fill(edge_tables)
    keys <- paste(edges$source, edges$target, edges$backend_database,
                  edges$statement_id, sep = "\r")
    groups <- split(seq_len(nrow(edges)), factor(keys, levels = unique(keys)))
    kept_rows <- vapply(groups, function(rows) {
        evidence_count <- edges$evidence_count[rows]
        evidence_count[is.na(evidence_count)] <- -Inf
        rows[which.max(evidence_count)]
    }, integer(1))
    query_types <- vapply(groups, function(rows) {
        .join_unique_values(edges$query_type[rows])
    }, character(1))
    copies_differ <- vapply(groups, function(rows) {
        length(unique(edges$evidence_count[rows])) > 1 ||
            length(unique(edges$confidence[rows])) > 1
    }, logical(1))
    if (any(copies_differ)) {
        message("Kept the copy with the most evidence for ",
                sum(copies_differ), " edge(s) whose copies differ in ",
                "evidence_count or confidence, statement_id: ",
                .list_values_for_message(
                    edges$statement_id[kept_rows[copies_differ]]), ".")
    }
    edges <- edges[kept_rows, , drop = FALSE]
    edges$query_type <- unname(query_types)
    rownames(edges) <- NULL
    edges
}

#' Combine nodes and collapse the rows of each node
#' @param node_tables list of nodes data.frames
#' @param entities \code{NULL}, or the entity table to recompute node
#' status from
#' @return nodes data.frame, one row per id and site
#' @keywords internal
#' @noRd
.merge_nodes <- function(node_tables, entities) {
    nodes <- .bind_rows_with_fill(node_tables)
    keys <- .build_node_keys(nodes$id, .get_column_or_default(nodes, "site",
                                                              NA_character_))
    groups <- split(seq_len(nrow(nodes)), factor(keys, levels = unique(keys)))
    groups <- groups[lengths(groups) > 1]
    if (is.null(entities)) {
        .check_node_status_agrees(nodes, groups)
    }
    merged <- nodes[!duplicated(keys), , drop = FALSE]
    repeated <- match(names(groups), keys[!duplicated(keys)])
    for (column in colnames(nodes)) {
        collapse <- switch(column,
                           node_role          = .join_unique_values,
                           has_measured_sites = .any_or_na,
                           included_in_query  = .any_or_na,
                           .first_non_na)
        values <- nodes[[column]]
        merged[[column]][repeated] <- vapply(groups, function(rows) {
            collapse(values[rows])
        }, values[NA_integer_], USE.NAMES = FALSE)
    }
    if (is.null(entities)) {
        merged <- .spread_measured_sites(merged)
    } else {
        merged <- .apply_entity_status(merged, entities)
    }
    rownames(merged) <- NULL
    merged
}

#' Give every row of a protein has_measured_sites = TRUE when any row has it
#'
#' has_measured_sites describes the protein, not the row. Merging a protein
#' network with a PTM network gives the protein a row from each, with
#' different sites, so they aren't collapsed into one row and can disagree.
#' @param nodes merged nodes data.frame
#' @return \code{nodes}
#' @keywords internal
#' @noRd
.spread_measured_sites <- function(nodes) {
    if (!"has_measured_sites" %in% colnames(nodes)) {
        return(nodes)
    }
    ids_with_sites <- nodes$id[nodes$has_measured_sites %in% TRUE]
    nodes$has_measured_sites[nodes$id %in% ids_with_sites] <- TRUE
    nodes
}

#' Stop when the copies of a node disagree on its status or statistics
#' @param nodes combined nodes data.frame
#' @param groups for each node with several rows, its row numbers
#' @return \code{NULL}, invisibly; errors otherwise
#' @keywords internal
#' @noRd
.check_node_status_agrees <- function(nodes, groups) {
    columns <- intersect(NODE_STATUS_COLUMNS, colnames(nodes))
    # One row per status column, one column per node with several rows
    column_differs <- matrix(vapply(groups, function(rows) {
        vapply(columns, function(column) {
            length(unique(nodes[[column]][rows])) > 1
        }, logical(1))
    }, logical(length(columns))), nrow = length(columns))
    disagrees <- colSums(column_differs) > 0
    if (any(disagrees)) {
        differing_columns <- columns[rowSums(column_differs) > 0]
        node_ids <- unique(nodes$id[vapply(groups[disagrees], `[`,
                                           integer(1), 1)])
        stop("The networks disagree on ",
             paste(differing_columns, collapse = ", "), " for node(s) ",
             .list_values_for_message(node_ids), ", so they were built ",
             "from different entity tables. Pass the entity ",
             "table as `entities =` to recompute node status from it.",
             call. = FALSE)
    }
    invisible(NULL)
}

#' Recompute measured, logFC, adj.pvalue, and has_measured_sites from an
#' entity table
#'
#' A node is measured when an entity row has its id (the parent protein's
#' id for PTM sites, see \code{.get_node_ids()}) and site. It has measured
#' sites when the entity table has a PTM site row for its id, the rule
#' get_network() uses.
#' @param nodes merged nodes data.frame
#' @param entities entity table
#' @return \code{nodes}, with \code{measured}, \code{logFC},
#' \code{adj.pvalue}, and \code{has_measured_sites} from \code{entities}
#' @keywords internal
#' @noRd
.apply_entity_status <- function(nodes, entities) {
    entity_keys <- .build_node_keys(.get_node_ids(entities),
                                    .get_column_or_default(entities, "site",
                                                           NA_character_))
    entity_rows <- match(.build_node_keys(
        nodes$id, .get_column_or_default(nodes, "site", NA_character_)),
        entity_keys)
    nodes$measured <- !is.na(entity_rows)
    for (column in c("logFC", "adj.pvalue")) {
        nodes[[column]] <- as.numeric(.get_column_or_default(
            entities, column, NA_real_)[entity_rows])
    }
    nodes$has_measured_sites <- nodes$id %in%
        .get_node_ids(entities)[entities$entity_type == "ptm_site"]
    nodes
}

#' One key per node: its id and site
#' @param ids node ids
#' @param sites sites, \code{NA} for nodes that aren't PTM sites
#' @return character vector
#' @keywords internal
#' @noRd
.build_node_keys <- function(ids, sites) {
    paste(ids, ifelse(is.na(sites), "", sites), sep = "\r")
}

#' Join the distinct values of a ;-joined column
#'
#' Splits each value on \code{";"}, so that already-merged values combine
#' too, and keeps the order of first appearance.
#' @param values character vector
#' @return a single string, or \code{NA} when every value is \code{NA}
#' @keywords internal
#' @noRd
.join_unique_values <- function(values) {
    values <- unlist(strsplit(values[!is.na(values)], ";", fixed = TRUE))
    if (length(values) == 0) {
        return(NA_character_)
    }
    paste(unique(values), collapse = ";")
}

#' TRUE if any value is TRUE, NA if every value is NA, otherwise FALSE
#' @keywords internal
#' @noRd
.any_or_na <- function(values) {
    if (all(is.na(values))) NA else any(values, na.rm = TRUE)
}

#' The first value that isn't NA, or NA
#' @keywords internal
#' @noRd
.first_non_na <- function(values) {
    present <- values[!is.na(values)]
    if (length(present) == 0) values[NA_integer_] else present[1]
}

#' Combine data.frames by row, filling missing columns with NA
#'
#' A column missing from a table is filled with \code{NA} of the type it
#' has in the first table that has it. Columns are in order of first
#' appearance.
#' @param tables list of data.frames
#' @return data.frame
#' @keywords internal
#' @noRd
.bind_rows_with_fill <- function(tables) {
    columns <- unique(unlist(lapply(tables, colnames)))
    templates <- lapply(columns, function(column) {
        for (table in tables) {
            if (column %in% colnames(table)) return(table[[column]])
        }
    })
    names(templates) <- columns
    filled <- lapply(tables, function(table) {
        table <- as.data.frame(table, stringsAsFactors = FALSE)
        for (column in setdiff(columns, colnames(table))) {
            table[[column]] <- templates[[column]][rep(NA_integer_,
                                                       nrow(table))]
        }
        table[, columns, drop = FALSE]
    })
    do.call(rbind, filled)
}
