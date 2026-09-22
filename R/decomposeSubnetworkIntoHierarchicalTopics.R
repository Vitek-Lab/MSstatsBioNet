#' Recursively decompose a subnetwork into a hierarchy of topic subnetworks
#'
#' Repeatedly applies \code{\link{decomposeSubnetworkByTopic}} to its own
#' topic subnetworks until every branch is small enough to inspect by hand
#' (at most \code{max_edges} edges), producing a topic tree: broad themes near
#' the root and increasingly specific sub-themes towards the leaves.
#'
#' INDRA evidence and PubMed abstracts are gathered once for the input
#' subnetwork and reused for every sub-decomposition, so no further network
#' requests are made during the recursion. Each sub-decomposition rebuilds its
#' vocabulary and refits the NMF on only the papers supporting that branch's
#' edges, which lets finer topics emerge.
#'
#' A branch stops splitting (becomes a leaf) when any of the following holds,
#' recorded in the \code{stop_reason} column of \code{tree}:
#' \describe{
#'   \item{small_enough}{it has at most \code{max_edges} edges.}
#'   \item{max_depth}{it sits at depth \code{max_depth}.}
#'   \item{too_few_papers}{fewer than two papers support its edges.}
#'   \item{no_split}{every child topic contained all of its edges (or none),
#'     so splitting would make no progress.}
#'   \item{failed: <message>}{the decomposition raised an error, e.g. no
#'     usable words in the branch's abstracts.}
#' }
#' Edges without any PMID-backed evidence cannot be assigned to a topic and
#' only appear at the root.
#'
#' @param subnetwork list with \code{nodes} and \code{edges} data.frames, e.g.
#'   the output of \code{\link{getSubnetworkFromIndra}}.
#' @param max_edges a branch with at most this many edges is not split
#'   further. Default 10.
#' @param n_topics number of topics per split. Default 5.
#' @param edge_topic_cutoff topic-share threshold for assigning an edge to a
#'   topic at each split; see \code{\link{decomposeSubnetworkByTopic}}. The
#'   default of 0.9 gives a near-partition, so sibling topics rarely share
#'   edges. Lower values allow an edge to appear in several sibling topics.
#' @param max_depth maximum depth of the tree (the root is depth 0). Default 5.
#' @param evidence optional pre-fetched evidence data.frame, e.g.
#'   \code{attr(topics, "corpus")$evidence} from
#'   \code{\link{decomposeSubnetworkByTopic}} or \code{result$corpus$evidence}
#'   from a previous call of this function. Default \code{NULL} queries INDRA
#'   once.
#' @param abstracts optional named character vector mapping PMID to abstract
#'   text. Only missing PMIDs are fetched from PubMed. Default \code{NULL}.
#' @param ... further arguments passed to
#'   \code{\link{decomposeSubnetworkByTopic}}, e.g. \code{n_top_terms},
#'   \code{min_term_count}, \code{include_ppi}, \code{seed}.
#'
#' @return An object of class \code{topicHierarchy}: a list with
#'   \describe{
#'     \item{tree}{data.frame with one row per topic, in depth-first order.
#'       Columns: \code{id} (\code{"root"}, \code{"1"}, \code{"1.2"}, ...),
#'       \code{parent_id}, \code{depth}, \code{topic} (index within the
#'       parent's split), \code{n_edges}, \code{n_nodes}, \code{n_papers}
#'       (papers supporting the topic's edges), \code{n_children},
#'       \code{is_leaf}, \code{stop_reason}, \code{mean_topic_weight} (mean
#'       share of the topic's edges' loading, a cohesion score),
#'       \code{top_terms} (collapsed with \code{", "}), \code{label}, and
#'       \code{pathString} (\code{"root/1/1.2"}, for \code{data.tree}).}
#'     \item{subnetworks}{named list keyed by \code{id}; each element is a
#'       subnetwork (\code{nodes}, \code{edges}, \code{topTerms},
#'       \code{pmids}) that can be passed to \code{\link{cytoscapeNetwork}} or
#'       \code{\link{exportNetworkToHTML}}.}
#'     \item{edge_membership}{long data.frame with one row per (topic, edge):
#'       \code{id}, \code{depth}, \code{is_leaf}, \code{source},
#'       \code{target}, \code{interaction}, \code{topicWeight}. Filter on
#'       \code{is_leaf} to see which fine-grained topic(s) each edge ends up
#'       in.}
#'     \item{corpus}{list with the \code{evidence} and \code{abstracts} used,
#'       for reuse via the \code{evidence} and \code{abstracts} arguments.}
#'     \item{params}{the settings used.}
#'   }
#'
#' @seealso \code{\link{decomposeSubnetworkByTopic}}
#' @note \strong{Beta feature:} This function is experimental and the API may
#'   change without notice in future versions.
#' @export
#'
#' @examples
#' \dontrun{
#' input <- data.table::fread(system.file(
#'     "extdata/groupComparisonModel.csv",
#'     package = "MSstatsBioNet"
#' ))
#' subnetwork <- getSubnetworkFromIndra(input)
#' hierarchy <- decomposeSubnetworkIntoHierarchicalTopics(
#'     subnetwork, max_edges = 10, n_topics = 5, edge_topic_cutoff = 0.9
#' )
#' hierarchy                       # indented topic tree
#' leaves <- hierarchy$tree[hierarchy$tree$is_leaf, ]
#' leaves[order(-leaves$mean_topic_weight), c("id", "n_edges", "top_terms")]
#'
#' # Inspect one fine-grained topic as a network.
#' leaf <- hierarchy$subnetworks[["1.2"]]
#' exportNetworkToHTML(leaf$nodes, leaf$edges)
#'
#' # Tree visualization with other packages, e.g.
#' # data.tree::as.Node(hierarchy$tree)
#' # igraph::graph_from_data_frame(
#' #     hierarchy$tree[-1, c("parent_id", "id")])
#' }
decomposeSubnetworkIntoHierarchicalTopics <- function(subnetwork,
                                                      max_edges = 10,
                                                      n_topics = 5,
                                                      edge_topic_cutoff = 0.9,
                                                      max_depth = 5,
                                                      evidence = NULL,
                                                      abstracts = NULL,
                                                      ...) {

    .validateDecomposeSubnetworkByTopicInput(subnetwork, n_topics,
                                             edge_topic_cutoff)
    .validateTopicCorpusInput(evidence, abstracts)
    if (!is.numeric(max_edges) || length(max_edges) != 1L ||
        is.na(max_edges) || max_edges < 1) {
        stop("`max_edges` must be a single number >= 1.")
    }
    if (!is.numeric(max_depth) || length(max_depth) != 1L ||
        is.na(max_depth) || max_depth < 0 ||
        max_depth != as.integer(max_depth)) {
        stop("`max_depth` must be a single non-negative integer.")
    }

    # Gather the corpus once; every sub-decomposition reuses it.
    corpus <- .gatherTopicCorpus(subnetwork$edges, evidence, abstracts)

    root <- list(nodes = subnetwork$nodes, edges = subnetwork$edges,
                 topTerms = character(0), pmids = character(0))
    stack <- list(list(id = "root", parent_id = NA_character_, depth = 0L,
                       topic = NA_integer_, subnetwork = root))
    rows <- list()
    subnetworks <- list()

    # Depth-first traversal so rows come out in tree (print) order.
    while (length(stack) > 0) {
        current <- stack[[length(stack)]]
        stack[[length(stack)]] <- NULL
        sub <- current$subnetwork

        split <- .splitTopicNode(sub, current$depth, corpus, max_edges,
                                 max_depth, n_topics, edge_topic_cutoff, ...)
        children <- split$children
        child_ids <- paste0(if (current$id == "root") "" else
                                paste0(current$id, "."),
                            names(children))
        for (i in rev(seq_along(children))) {
            stack[[length(stack) + 1]] <- list(
                id = child_ids[i], parent_id = current$id,
                depth = current$depth + 1L, topic = children[[i]]$topic,
                subnetwork = children[[i]]
            )
        }

        rows[[length(rows) + 1]] <- .topicNodeRow(current, split,
                                                  length(children))
        subnetworks[[current$id]] <- sub[c("nodes", "edges",
                                           "topTerms", "pmids")]
    }

    tree <- do.call(rbind, rows)
    rownames(tree) <- NULL
    tree$pathString <- .topicPathStrings(tree$id, tree$parent_id)

    result <- list(
        tree            = tree,
        subnetworks     = subnetworks,
        edge_membership = .topicEdgeMembership(tree, subnetworks),
        corpus          = corpus,
        params          = list(max_edges = max_edges, n_topics = n_topics,
                               edge_topic_cutoff = edge_topic_cutoff,
                               max_depth = max_depth, ...)
    )
    class(result) <- "topicHierarchy"
    result
}


#' Print a topic hierarchy as an indented tree
#'
#' @param x a \code{topicHierarchy} object from
#'   \code{\link{decomposeSubnetworkIntoHierarchicalTopics}}.
#' @param n_terms number of top terms to show per topic. Default 5.
#' @param ... ignored.
#' @return \code{x}, invisibly.
#' @export
print.topicHierarchy <- function(x, n_terms = 5, ...) {
    tree <- x$tree
    cat(sprintf(
        "Topic hierarchy: %d topics, %d leaves, depth %d (max_edges = %s)\n",
        nrow(tree) - 1L, sum(tree$is_leaf & tree$depth > 0),
        max(tree$depth), format(x$params$max_edges)
    ))
    for (i in seq_len(nrow(tree))) {
        terms <- x$subnetworks[[tree$id[i]]]$topTerms
        terms <- paste(utils::head(terms, n_terms), collapse = ", ")
        stop_note <- if (tree$is_leaf[i] && tree$depth[i] > 0 &&
                         tree$stop_reason[i] != "small_enough") {
            paste0("  <", tree$stop_reason[i], ">")
        } else ""
        cat(sprintf("%s%s [%d edges, %d papers]%s%s\n",
                    strrep("  ", tree$depth[i]), tree$id[i],
                    tree$n_edges[i], tree$n_papers[i],
                    if (nzchar(terms)) paste0(" ", terms) else "",
                    stop_note))
    }
    invisible(x)
}


#' Decide whether a topic node should split, and split it if so
#'
#' @param sub topic subnetwork (list with `nodes` and `edges`)
#' @param depth depth of the node in the tree
#' @param corpus list with `evidence` and `abstracts` from
#'   \code{.gatherTopicCorpus}
#' @param max_edges,max_depth,n_topics,edge_topic_cutoff see
#'   \code{decomposeSubnetworkIntoHierarchicalTopics}
#' @param ... passed to \code{decomposeSubnetworkByTopic}
#' @return list with `children` (named list of child topic subnetworks, empty
#'   for a leaf), `stop_reason` (NA when the node split), and `n_papers`
#' @keywords internal
#' @noRd
.splitTopicNode <- function(sub, depth, corpus, max_edges, max_depth,
                            n_topics, edge_topic_cutoff, ...) {
    n_edges <- nrow(sub$edges)
    node_evidence <- .subsetEvidenceToEdges(corpus$evidence, sub$edges)
    n_papers <- length(unique(node_evidence$pmid))
    leaf <- function(reason) {
        list(children = list(), stop_reason = reason, n_papers = n_papers)
    }

    if (n_edges <= max_edges) return(leaf("small_enough"))
    if (depth >= max_depth) return(leaf("max_depth"))
    if (n_papers < 2) return(leaf("too_few_papers"))

    topics <- tryCatch(
        withCallingHandlers(
            decomposeSubnetworkByTopic(
                sub, n_topics = n_topics,
                edge_topic_cutoff = edge_topic_cutoff,
                evidence = node_evidence, abstracts = corpus$abstracts, ...
            ),
            # n_topics is routinely capped by the paper count deep in the tree.
            warning = function(w) {
                if (grepl("^Only \\d+ papers available",
                          conditionMessage(w))) {
                    invokeRestart("muffleWarning")
                }
            }
        ),
        error = function(e) e
    )
    if (inherits(topics, "error")) {
        return(leaf(paste0("failed: ", conditionMessage(topics))))
    }

    # Keep only children that strictly shrink the branch; this also
    # guarantees the recursion terminates.
    n_child_edges <- vapply(topics, function(t) nrow(t$edges), integer(1))
    children <- topics[n_child_edges > 0 & n_child_edges < n_edges]
    if (length(children) == 0) return(leaf("no_split"))
    names(children) <- vapply(children, function(t) as.character(t$topic),
                              character(1))

    list(children = children, stop_reason = NA_character_,
         n_papers = n_papers)
}


#' Build the summary row of the tree data.frame for one topic node
#' @param current stack entry (id, parent_id, depth, topic, subnetwork)
#' @param split output of \code{.splitTopicNode}
#' @param n_children number of children the node was split into
#' @return one-row data.frame
#' @keywords internal
#' @noRd
.topicNodeRow <- function(current, split, n_children) {
    sub <- current$subnetwork
    weights <- sub$edges$topicWeight
    top_terms <- paste(sub$topTerms, collapse = ", ")
    data.frame(
        id                = current$id,
        parent_id         = current$parent_id,
        depth             = current$depth,
        topic             = current$topic,
        n_edges           = nrow(sub$edges),
        n_nodes           = nrow(sub$nodes),
        n_papers          = split$n_papers,
        n_children        = n_children,
        is_leaf           = n_children == 0,
        stop_reason       = split$stop_reason,
        mean_topic_weight = if (is.null(weights) || current$id == "root")
                                NA_real_ else mean(weights),
        top_terms         = top_terms,
        label             = sprintf("%s (%d edges)%s", current$id,
                                    nrow(sub$edges),
                                    if (nzchar(top_terms))
                                        paste0(": ", paste(utils::head(
                                            sub$topTerms, 3),
                                            collapse = ", "))
                                    else ""),
        stringsAsFactors  = FALSE
    )
}


#' Build data.tree-style path strings ("root/1/1.2") for each tree node
#' @param ids character vector of node ids
#' @param parent_ids character vector of parent ids (NA for the root)
#' @return character vector of path strings aligned to `ids`
#' @keywords internal
#' @noRd
.topicPathStrings <- function(ids, parent_ids) {
    parent_of <- stats::setNames(parent_ids, ids)
    vapply(ids, function(id) {
        path <- id
        while (!is.na(parent_of[[id]])) {
            id <- parent_of[[id]]
            path <- c(id, path)
        }
        paste(path, collapse = "/")
    }, character(1), USE.NAMES = FALSE)
}


#' Long table of topic-to-edge membership across the hierarchy
#' @param tree tree data.frame
#' @param subnetworks named list of topic subnetworks keyed by tree id
#' @return data.frame with one row per (non-root topic, edge)
#' @keywords internal
#' @noRd
.topicEdgeMembership <- function(tree, subnetworks) {
    non_root <- tree[tree$id != "root", , drop = FALSE]
    parts <- lapply(seq_len(nrow(non_root)), function(i) {
        edges <- subnetworks[[non_root$id[i]]]$edges
        if (nrow(edges) == 0) return(NULL)
        data.frame(
            id          = non_root$id[i],
            depth       = non_root$depth[i],
            is_leaf     = non_root$is_leaf[i],
            source      = edges$source,
            target      = edges$target,
            interaction = edges$interaction,
            topicWeight = edges$topicWeight,
            stringsAsFactors = FALSE
        )
    })
    parts <- parts[!vapply(parts, is.null, logical(1))]
    if (length(parts) == 0) {
        return(data.frame(id = character(), depth = integer(),
                          is_leaf = logical(), source = character(),
                          target = character(), interaction = character(),
                          topicWeight = numeric(),
                          stringsAsFactors = FALSE))
    }
    out <- do.call(rbind, parts)
    rownames(out) <- NULL
    out
}
