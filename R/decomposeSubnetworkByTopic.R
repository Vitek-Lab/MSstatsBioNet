#' Decompose a subnetwork into topic-specific subnetworks via joint NMF
#'
#' Takes a subnetwork (the output of \code{\link{getSubnetworkFromIndra}}) and
#' splits it into a list of smaller, topic-specific subnetworks discovered with
#' unsupervised non-negative matrix factorization (NMF).
#'
#' The procedure is:
#' \enumerate{
#'   \item For every edge, the supporting INDRA evidence is retrieved and the
#'     PubMed abstract of each referenced PMID is fetched. Papers (PMIDs) are
#'     the shared unit of analysis.
#'   \item Two matrices are built that share the same rows (papers):
#'     \code{X_text} (papers x words, term counts from the abstracts) and
#'     \code{X_edges} (papers x unique \code{source_target_interaction}
#'     combinations, evidence-sentence counts).
#'   \item NMF learns a basis matrix \code{W} (papers x topics). When
#'     \code{include_ppi = TRUE} (the default) a \emph{joint} NMF learns a
#'     single shared \code{W} such that \eqn{X_{text} \approx W H_{text}} and
#'     \eqn{X_{edges} \approx W H_{edges}}, tying each learned topic to both a
#'     set of words and a set of edges. When \code{include_ppi = FALSE} the
#'     factorization uses only \code{X_text} (\eqn{X_{text} \approx W H_{text}});
#'     the PPI evidence is excluded from the modeling and edge-topic loadings
#'     are instead derived afterwards by folding the edge counts onto the
#'     text-learned topics (\eqn{H_{edges} = W^\top X_{edges}}). This lets you
#'     compare topic structure with and without the PPI view.
#'   \item Each topic becomes its own subnetwork: an edge is included in a
#'     topic when that topic carries at least \code{edge_topic_cutoff} of the
#'     edge's loading (soft, overlapping assignment), and nodes are restricted
#'     to those touched by the kept edges.
#' }
#'
#' @param subnetwork list with \code{nodes} and \code{edges} data.frames, e.g.
#'   the output of \code{\link{getSubnetworkFromIndra}}.
#' @param n_topics number of topics (rank of the factorization). Default 5.
#' @param edge_topic_cutoff numeric in \code{[0, 1]}; an edge is added to a
#'   topic's subnetwork when the topic carries at least this share of the
#'   edge's total loading. Each edge is always included in at least its
#'   highest-loading topic. Default 0.2.
#' @param n_top_terms number of top words to report per topic. Default 10.
#' @param min_term_count minimum corpus frequency for a word to be kept when
#'   building \code{X_text}. Default 2.
#' @param max_iter maximum number of NMF multiplicative-update iterations.
#'   Default 200.
#' @param tol relative-change tolerance for NMF early stopping. Default 1e-4.
#' @param seed random seed for NMF initialization. Default 1.
#' @param include_ppi logical; if \code{TRUE} (default) the PPI/edge matrix is
#'   factorized jointly with the text matrix via a shared basis. If
#'   \code{FALSE}, NMF is run on the paper-word matrix only and edge-topic
#'   loadings are derived afterwards by folding edge counts onto the
#'   text-learned topics, so the PPIs do not influence the topics themselves.
#'
#' @return A list of length \code{n_topics}, named \code{topic_1} ...
#'   \code{topic_k}. Each element is a topic-specific subnetwork: a list with
#'   \describe{
#'     \item{nodes}{nodes data.frame restricted to the topic's edges.}
#'     \item{edges}{edges data.frame for the topic, with an added
#'       \code{topicWeight} column (the edge's topic share).}
#'     \item{topic}{the topic index.}
#'     \item{topTerms}{character vector of the topic's top words.}
#'     \item{pmids}{PMIDs whose strongest topic loading is this topic.}
#'   }
#'   The full factorization (W, H_text, H_edges, etc.) is attached as the
#'   \code{"nmf"} attribute of the returned list.
#'
#' @seealso \code{\link{getSubnetworkFromIndra}},
#'   \code{\link{filterSubnetworkByContext}}
#'
#' @export
#' 
#' @note **Beta feature:** This function is experimental and the API may
#'   change without notice in future versions.
#'
#' @examples
#' \dontrun{
#' input <- data.table::fread(system.file(
#'     "extdata/groupComparisonModel.csv",
#'     package = "MSstatsBioNet"
#' ))
#' subnetwork <- getSubnetworkFromIndra(input)
#' topics <- decomposeSubnetworkByTopic(subnetwork, n_topics = 5)
#' topics$topic_1$topTerms
#' exportNetworkToHTML(topics$topic_1$nodes, topics$topic_1$edges)
#' }
decomposeSubnetworkByTopic <- function(subnetwork,
                                       n_topics = 5,
                                       edge_topic_cutoff = 0.2,
                                       n_top_terms = 10,
                                       min_term_count = 2,
                                       max_iter = 200,
                                       tol = 1e-4,
                                       seed = 1,
                                       include_ppi = TRUE) {

    .validateDecomposeSubnetworkByTopicInput(subnetwork, n_topics,
                                             edge_topic_cutoff, include_ppi)

    # 1-3. Build the shared paper-by-word and paper-by-edge matrices.
    mats <- .buildTopicMatrices(subnetwork, n_topics, min_term_count)
    nodes     <- mats$nodes
    edges     <- mats$edges
    pmids     <- mats$pmids
    edge_keys <- mats$edge_keys
    X_text    <- mats$X_text
    X_edges   <- mats$X_edges
    n_topics  <- mats$n_topics

    # 4. NMF: jointly over text + PPIs, or over text only.
    if (include_ppi) {
        model <- .jointNMF(X_text, X_edges, k = n_topics,
                           max_iter = max_iter, tol = tol, seed = seed)
    } else {
        model <- .textNMF(X_text, k = n_topics,
                          max_iter = max_iter, tol = tol, seed = seed)
        # Edges are excluded from the modeling, but still need topic loadings:
        # fold the edge counts onto the text-learned topics.
        model$H_edges <- .edgeLoadingsFromTopics(model$W, X_edges)
    }

    # 5. One subnetwork per topic (soft / overlapping edge assignment).
    shares <- .edgeTopicShares(model$H_edges)            # topics x edges
    edge_argmax <- apply(shares, 2, which.max)           # best topic per edge
    edges_key_vec <- .edgeKey(edges$source, edges$target, edges$interaction)
    paper_argmax <- apply(model$W, 1, which.max)         # best topic per paper

    topics <- lapply(seq_len(n_topics), function(t) {
        # Edge keys assigned to this topic.
        topic_keys <- edge_keys[shares[t, ] >= edge_topic_cutoff |
                                    edge_argmax == t]
        in_topic <- edges_key_vec %in% topic_keys
        topic_edges <- edges[in_topic, , drop = FALSE]
        if (nrow(topic_edges) > 0) {
            topic_edges$topicWeight <-
                shares[t, match(edges_key_vec[in_topic], edge_keys)]
        }
        topic_nodes <- nodes[nodes$id %in%
                                 c(topic_edges$source, topic_edges$target), ,
                             drop = FALSE]
        list(
            nodes    = topic_nodes,
            edges    = topic_edges,
            topic    = t,
            topTerms = .topTermsForTopic(model$H_text, t, n_top_terms),
            pmids    = pmids[paper_argmax == t]
        )
    })
    names(topics) <- paste0("topic_", seq_len(n_topics))

    attr(topics, "nmf") <- list(
        W           = model$W,
        H_text      = model$H_text,
        H_edges     = model$H_edges,
        terms       = colnames(X_text),
        edge_keys   = edge_keys,
        pmids       = pmids,
        objective   = model$objective,
        n_iter      = model$n_iter,
        include_ppi = include_ppi
    )
    return(topics)
}
