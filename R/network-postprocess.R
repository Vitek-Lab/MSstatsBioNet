# Filters that run on a finished network, after get_network().

#' Remove evidence curated as incorrect
#'
#' Subtracts the evidence curated as incorrect, from
#' \code{\link{get_curations}()}, from each edge's \code{evidence_count},
#' then drops the edges left below \code{min_evidence}, and the nodes those
#' edges leave without any edge. Edges left with no evidence are dropped
#' whatever \code{min_evidence} is. A message says how many edges were
#' dropped. Edges with no evidence count (\code{NA}), such as STRING's,
#' are kept as they are, and their backend is not asked for curations.
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
    # Edges with no evidence count (NA, e.g. STRING) have nothing to subtract
    # from, so they are kept and their backend is not asked
    counted <- !is.na(edges$evidence_count)
    incorrect_counts <- .count_incorrect_evidence(
        edges[counted, , drop = FALSE], backend)
    edges$evidence_count[counted] <- as.integer(
        edges$evidence_count[counted] - incorrect_counts)
    edges$evidence_count <- as.integer(edges$evidence_count)
    keep <- !counted |
        (edges$evidence_count >= min_evidence & edges$evidence_count >= 1)
    if (any(!keep)) {
        message("Dropping ", sum(!keep), " edge(s) with fewer than ",
                max(min_evidence, 1), " evidence after removing the ",
                "evidence curated as incorrect.")
    }
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
