#' Validate input for decomposeSubnetworkByTopic
#' @param subnetwork list with `nodes` and `edges` data.frames
#' @param n_topics number of topics (rank of the factorization)
#' @param edge_topic_cutoff topic-share threshold for assigning an edge to a topic
#' @keywords internal
#' @noRd
.validateDecomposeSubnetworkByTopicInput <- function(subnetwork,
                                                     n_topics,
                                                     edge_topic_cutoff) {
    if (!is.list(subnetwork) ||
        !all(c("nodes", "edges") %in% names(subnetwork))) {
        stop("`subnetwork` must be a list containing `nodes` and `edges`, ",
             "e.g. the output of getSubnetworkFromIndra().")
    }
    if (!is.data.frame(subnetwork$nodes) || !"id" %in% names(subnetwork$nodes)) {
        stop("`subnetwork$nodes` must be a data.frame with an `id` column.")
    }
    required_edge_cols <- c("source", "target", "interaction",
                            "site", "evidenceLink", "stmt_hash")
    missing_cols <- setdiff(required_edge_cols, names(subnetwork$edges))
    if (!is.data.frame(subnetwork$edges) || length(missing_cols) > 0) {
        stop(sprintf(
            "`subnetwork$edges` must be a data.frame with columns: %s",
            paste(required_edge_cols, collapse = ", ")
        ))
    }
    if (!is.numeric(n_topics) || length(n_topics) != 1L || is.na(n_topics) ||
        n_topics < 1 || n_topics != as.integer(n_topics)) {
        stop("`n_topics` must be a single positive integer.")
    }
    if (!is.numeric(edge_topic_cutoff) || length(edge_topic_cutoff) != 1L ||
        is.na(edge_topic_cutoff) || edge_topic_cutoff < 0 ||
        edge_topic_cutoff > 1) {
        stop("`edge_topic_cutoff` must be a single numeric value in [0, 1].")
    }
}


#' Build a unique source_target_interaction key for each edge
#' @param source character vector of source node ids
#' @param target character vector of target node ids
#' @param interaction character vector of interaction types
#' @return character vector of edge keys
#' @keywords internal
#' @noRd
.edgeKey <- function(source, target, interaction) {
    paste(source, target, interaction, sep = "||")
}


#' Build the paper-by-word matrix (X_text) from PubMed abstracts
#'
#' Tokenises abstracts with text2vec and returns a dense paper-by-word count
#' matrix whose rows are aligned to `pmids`.
#'
#' @param pmids character vector of PubMed IDs (rows of the matrix)
#' @param abstracts character vector of abstract texts, aligned to `pmids`
#' @param min_term_count minimum corpus term frequency to keep a word
#' @return dense numeric matrix (papers x words) with rownames = pmids
#' @keywords internal
#' @noRd
#' @importFrom text2vec itoken word_tokenizer create_vocabulary
#'   prune_vocabulary vocab_vectorizer create_dtm
#' @importFrom stopwords stopwords
.buildTextMatrix <- function(pmids, abstracts, min_term_count = 2) {
    tokens <- itoken(abstracts,
                     preprocessor = tolower,
                     tokenizer    = word_tokenizer,
                     ids          = pmids,
                     progressbar  = FALSE)
    vocab <- create_vocabulary(tokens, stopwords = stopwords("en"))
    pruned <- prune_vocabulary(vocab, term_count_min = min_term_count)
    if (nrow(pruned) == 0) {
        # Fall back to keeping every term before giving up.
        pruned <- prune_vocabulary(vocab, term_count_min = 1)
    }
    if (nrow(pruned) == 0) {
        stop("No usable words found in the fetched abstracts; ",
             "cannot build the text matrix.")
    }
    vectorizer <- vocab_vectorizer(pruned)
    dtm <- create_dtm(tokens, vectorizer)
    dtm <- as.matrix(dtm)
    # Align rows to the requested pmid order (zero rows for empty abstracts).
    aligned <- matrix(0, nrow = length(pmids), ncol = ncol(dtm),
                      dimnames = list(pmids, colnames(dtm)))
    common <- intersect(pmids, rownames(dtm))
    if (length(common) > 0) {
        aligned[common, ] <- dtm[common, , drop = FALSE]
    }
    return(aligned)
}


#' Build the paper-by-edge matrix (X_edges) of evidence counts
#'
#' @param evidence data.frame from \code{.extract_evidence_text}; one row per
#'   evidence sentence with `pmid`, `source`, `target`, `interaction`.
#' @param pmids character vector of PubMed IDs (rows of the matrix)
#' @param edge_keys character vector of unique edge keys (columns of the matrix)
#' @return dense numeric matrix (papers x edges) of evidence-sentence counts
#' @keywords internal
#' @noRd
.buildEdgeMatrix <- function(evidence, pmids, edge_keys) {
    X <- matrix(0, nrow = length(pmids), ncol = length(edge_keys),
                dimnames = list(pmids, edge_keys))
    ev_key <- .edgeKey(evidence$source, evidence$target, evidence$interaction)
    counts <- table(factor(evidence$pmid, levels = pmids),
                    factor(ev_key, levels = edge_keys))
    X[] <- as.numeric(counts)
    return(X)
}


#' Joint non-negative matrix factorization with a shared basis matrix
#'
#' Factorizes two matrices that share the same rows (papers) so that
#' \eqn{X_{text} \approx W H_{text}} and \eqn{X_{edges} \approx W H_{edges}},
#' where the basis matrix \eqn{W} (papers x topics) is shared between both
#' views. Uses multiplicative updates that minimise the combined squared
#' Frobenius reconstruction error.
#'
#' @param X_text paper-by-word matrix
#' @param X_edges paper-by-edge matrix
#' @param k number of topics (rank)
#' @param max_iter maximum number of multiplicative-update iterations
#' @param tol relative-change tolerance for early stopping
#' @param seed random seed for initialization
#' @param normalize logical; if TRUE each view is scaled to unit Frobenius norm
#'   so neither view dominates the shared factorization
#' @return list with elements W, H_text, H_edges, objective (per-iteration
#'   objective values), and n_iter
#' @keywords internal
#' @noRd
#' @importFrom stats runif
.jointNMF <- function(X_text, X_edges, k,
                      max_iter = 200, tol = 1e-4, seed = 1,
                      normalize = TRUE) {
    eps <- 1e-10
    n <- nrow(X_text)
    if (nrow(X_edges) != n) {
        stop("X_text and X_edges must have the same number of rows (papers).")
    }
    if (normalize) {
        ft <- sqrt(sum(X_text^2))
        fe <- sqrt(sum(X_edges^2))
        if (ft > 0) X_text <- X_text / ft
        if (fe > 0) X_edges <- X_edges / fe
    }

    p <- ncol(X_text)
    m <- ncol(X_edges)

    set.seed(seed)
    W <- matrix(runif(n * k), n, k)
    H_text <- matrix(runif(k * p), k, p)
    H_edges <- matrix(runif(k * m), k, m)

    objective <- numeric(0)
    prev_obj <- Inf
    n_iter <- 0L
    for (iter in seq_len(max_iter)) {
        n_iter <- iter
        # Update view-specific coefficient matrices.
        WtW <- crossprod(W)                      # k x k
        H_text  <- H_text  * (crossprod(W, X_text))  / (WtW %*% H_text  + eps)
        H_edges <- H_edges * (crossprod(W, X_edges)) / (WtW %*% H_edges + eps)

        # Update the shared basis from both views jointly.
        numer <- tcrossprod(X_text, H_text) + tcrossprod(X_edges, H_edges)
        denom <- W %*% (tcrossprod(H_text) + tcrossprod(H_edges)) + eps
        W <- W * numer / denom

        obj <- sum((X_text - W %*% H_text)^2) +
            sum((X_edges - W %*% H_edges)^2)
        objective <- c(objective, obj)
        if (is.finite(prev_obj) &&
            abs(prev_obj - obj) <= tol * (prev_obj + eps)) {
            break
        }
        prev_obj <- obj
    }

    rownames(W) <- rownames(X_text)
    colnames(W) <- paste0("topic_", seq_len(k))
    rownames(H_text) <- colnames(W)
    colnames(H_text) <- colnames(X_text)
    rownames(H_edges) <- colnames(W)
    colnames(H_edges) <- colnames(X_edges)

    return(list(W = W, H_text = H_text, H_edges = H_edges,
                objective = objective, n_iter = n_iter))
}


#' Compute the per-edge topic shares from H_edges
#'
#' Column-normalises H_edges so that, for each edge, the loadings across topics
#' sum to 1, giving each edge a distribution over topics.
#'
#' @param H_edges topics-by-edges coefficient matrix
#' @return topics-by-edges matrix of topic shares (columns sum to 1)
#' @keywords internal
#' @noRd
.edgeTopicShares <- function(H_edges) {
    col_sums <- colSums(H_edges)
    col_sums[col_sums == 0] <- 1
    sweep(H_edges, 2, col_sums, "/")
}


#' Top words for a topic from H_text
#' @param H_text topics-by-words coefficient matrix
#' @param topic integer topic index
#' @param n number of top terms to return
#' @return character vector of the top-weighted terms for the topic
#' @keywords internal
#' @noRd
.topTermsForTopic <- function(H_text, topic, n = 10) {
    weights <- H_text[topic, ]
    weights <- weights[weights > 0]
    if (length(weights) == 0) return(character(0))
    ordered <- sort(weights, decreasing = TRUE)
    names(ordered)[seq_len(min(n, length(ordered)))]
}
