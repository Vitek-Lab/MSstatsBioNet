#' Validate input for decomposeSubnetworkByTopic
#' @param subnetwork list with `nodes` and `edges` data.frames
#' @param n_topics number of topics (rank of the factorization)
#' @param edge_topic_cutoff topic-share threshold for assigning an edge to a topic
#' @param include_ppi logical; whether the PPI/edge matrix is included in the NMF
#' @keywords internal
#' @noRd
.validateDecomposeSubnetworkByTopicInput <- function(subnetwork,
                                                     n_topics,
                                                     edge_topic_cutoff,
                                                     include_ppi = TRUE) {
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
    if (!is.logical(include_ppi) || length(include_ppi) != 1L ||
        is.na(include_ppi)) {
        stop("`include_ppi` must be a single logical value (TRUE or FALSE).")
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


#' Build the shared paper-by-word and paper-by-edge matrices for topic modeling
#'
#' Performs the data-preparation steps shared by
#' \code{\link{decomposeSubnetworkByTopic}} and
#' \code{\link{compareTopicModels}}: extracts INDRA evidence, fetches and
#' tokenises the PubMed abstracts, and builds the aligned \code{X_text}
#' (papers x words) and \code{X_edges} (papers x edges) matrices. Factored out
#' so the (network-bound) abstract fetch and matrix construction happen once and
#' can be reused across many NMF fits.
#'
#' @param subnetwork list with `nodes` and `edges` data.frames
#' @param n_topics requested number of topics; reduced (with a warning) when
#'   fewer papers than topics are available
#' @param min_term_count minimum corpus term frequency to keep a word
#' @return list with `nodes`, `edges`, `evidence`, `pmids`, `edge_keys`,
#'   `X_text`, `X_edges`, and the (possibly reduced) `n_topics`
#' @keywords internal
#' @noRd
.buildTopicMatrices <- function(subnetwork, n_topics, min_term_count = 2) {
    nodes <- subnetwork$nodes
    edges <- subnetwork$edges
    n_topics <- as.integer(n_topics)

    # 1. Evidence (paper <-> edge links) for every edge.
    evidence <- .extract_evidence_text(edges)
    evidence <- evidence[!is.na(evidence$pmid) & nchar(evidence$pmid) > 0, ]
    if (nrow(evidence) == 0) {
        stop("No evidence with PMIDs was found for any edge; ",
             "cannot decompose into topics.")
    }

    pmids <- unique(evidence$pmid)
    edge_keys <- unique(.edgeKey(evidence$source, evidence$target,
                                 evidence$interaction))

    if (length(pmids) < n_topics) {
        warning(sprintf(
            "Only %d papers available; reducing n_topics from %d to %d.",
            length(pmids), n_topics, length(pmids)
        ))
        n_topics <- length(pmids)
    }

    # 2. X_text (papers x words) from PubMed abstracts.
    abstract_list <- .fetch_clean_abstracts_xml(pmids)
    abstracts <- vapply(pmids, function(p) {
        a <- abstract_list[[p]]
        if (is.null(a)) "" else a
    }, character(1))
    X_text <- .buildTextMatrix(pmids, abstracts, min_term_count)

    # 3. X_edges (papers x source_target_interaction) of evidence counts.
    X_edges <- .buildEdgeMatrix(evidence, pmids, edge_keys)

    list(nodes = nodes, edges = edges, evidence = evidence,
         pmids = pmids, edge_keys = edge_keys,
         X_text = X_text, X_edges = X_edges, n_topics = n_topics)
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
#' @importFrom text2vec itoken word_tokenizer create_vocabulary prune_vocabulary vocab_vectorizer create_dtm
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


#' Single-view non-negative matrix factorization of the paper-by-word matrix
#'
#' Factorizes the text view only, \eqn{X_{text} \approx W H_{text}}, using the
#' same multiplicative updates and normalization as \code{\link{.jointNMF}} but
#' without a PPI/edge view. Used by \code{decomposeSubnetworkByTopic} when
#' \code{include_ppi = FALSE} so the learned topics reflect paper words alone.
#'
#' @param X_text paper-by-word matrix
#' @param k number of topics (rank)
#' @param max_iter maximum number of multiplicative-update iterations
#' @param tol relative-change tolerance for early stopping
#' @param seed random seed for initialization
#' @param normalize logical; if TRUE the view is scaled to unit Frobenius norm
#' @return list with elements W, H_text, objective (per-iteration objective
#'   values), and n_iter
#' @keywords internal
#' @noRd
#' @importFrom stats runif
.textNMF <- function(X_text, k,
                     max_iter = 200, tol = 1e-4, seed = 1,
                     normalize = TRUE) {
    eps <- 1e-10
    n <- nrow(X_text)
    if (normalize) {
        ft <- sqrt(sum(X_text^2))
        if (ft > 0) X_text <- X_text / ft
    }

    p <- ncol(X_text)

    set.seed(seed)
    W <- matrix(runif(n * k), n, k)
    H_text <- matrix(runif(k * p), k, p)

    objective <- numeric(0)
    prev_obj <- Inf
    n_iter <- 0L
    for (iter in seq_len(max_iter)) {
        n_iter <- iter
        WtW <- crossprod(W)                          # k x k
        H_text <- H_text * (crossprod(W, X_text)) / (WtW %*% H_text + eps)

        numer <- tcrossprod(X_text, H_text)
        denom <- W %*% tcrossprod(H_text) + eps
        W <- W * numer / denom

        obj <- sum((X_text - W %*% H_text)^2)
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

    return(list(W = W, H_text = H_text,
                objective = objective, n_iter = n_iter))
}


#' Derive topic-by-edge loadings from paper topics (text-only model)
#'
#' When PPIs are excluded from the factorization there is no learned
#' \code{H_edges}. Edges still need topic loadings so they can be assigned to
#' topic-specific subnetworks, so each edge's loading on a topic is obtained by
#' folding its evidence counts onto the text-learned paper topics:
#' \eqn{H_{edges} = W^\top X_{edges}}. The resulting matrix matches the shape
#' and dimnames of a jointly-learned \code{H_edges} and is consumed downstream
#' identically (e.g. by \code{\link{.edgeTopicShares}}).
#'
#' @param W paper-by-topic basis matrix from \code{\link{.textNMF}}
#' @param X_edges paper-by-edge matrix of evidence counts
#' @return topics-by-edges matrix of (unnormalized) edge-topic loadings
#' @keywords internal
#' @noRd
.edgeLoadingsFromTopics <- function(W, X_edges) {
    H_edges <- crossprod(W, X_edges)                 # topics x edges
    rownames(H_edges) <- colnames(W)
    colnames(H_edges) <- colnames(X_edges)
    return(H_edges)
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


#' Extract a hard topic assignment (partition) from a fitted NMF model
#'
#' Reduces a factorization to one integer label per unit so that solutions from
#' different runs can be compared with permutation-invariant metrics. For
#' \code{unit = "papers"} each paper is assigned its strongest topic in
#' \code{W}; for \code{unit = "edges"} each edge is assigned its strongest topic
#' share in \code{H_edges} (matching the edge assignment used by
#' \code{decomposeSubnetworkByTopic}).
#'
#' @param model fitted model list with `W` and `H_edges`
#' @param unit either "edges" or "papers"
#' @return integer vector of topic labels, named by edge key or pmid
#' @keywords internal
#' @noRd
.topicPartition <- function(model, unit = c("edges", "papers")) {
    unit <- match.arg(unit)
    if (unit == "papers") {
        labels <- apply(model$W, 1, which.max)
        names(labels) <- rownames(model$W)
    } else {
        shares <- .edgeTopicShares(model$H_edges)
        labels <- apply(shares, 2, which.max)
        names(labels) <- colnames(model$H_edges)
    }
    as.integer(labels)
}


#' Adjusted Rand Index between two partitions
#'
#' Permutation-invariant agreement between two labelings of the same items,
#' corrected for chance. Equals 1 for identical partitions and ~0 for the
#' agreement expected at random. Implemented from the contingency table so no
#' extra package dependency is needed.
#'
#' @param a,b integer/factor vectors of equal length (cluster labels)
#' @return numeric scalar in roughly \code{[-0.5, 1]}
#' @keywords internal
#' @noRd
.adjustedRandIndex <- function(a, b) {
    if (length(a) != length(b)) {
        stop("`a` and `b` must have the same length.")
    }
    n <- length(a)
    if (n == 0) return(NA_real_)
    choose2 <- function(x) x * (x - 1) / 2
    tab <- table(a, b)
    sum_cells <- sum(choose2(tab))
    sum_a <- sum(choose2(rowSums(tab)))
    sum_b <- sum(choose2(colSums(tab)))
    expected <- sum_a * sum_b / choose2(n)
    max_index <- (sum_a + sum_b) / 2
    if (max_index == expected) return(1)        # both trivially one cluster
    (sum_cells - expected) / (max_index - expected)
}


#' Consensus (co-membership) matrix across a set of partitions
#'
#' Entry \code{[i, j]} is the fraction of partitions in which items \code{i} and
#' \code{j} were assigned to the same topic. All partitions must label the same
#' items in the same order.
#'
#' @param partitions list of integer label vectors of equal length
#' @return symmetric numeric matrix with entries in \code{[0, 1]}
#' @keywords internal
#' @noRd
.consensusMatrix <- function(partitions) {
    n <- length(partitions[[1]])
    consensus <- matrix(0, n, n,
                        dimnames = list(names(partitions[[1]]),
                                        names(partitions[[1]])))
    for (p in partitions) {
        consensus <- consensus + outer(p, p, "==")
    }
    consensus / length(partitions)
}


#' Dispersion coefficient of a consensus matrix
#'
#' Summarises how close a consensus matrix is to perfectly stable (every entry
#' 0 or 1) versus maximally ambiguous (every entry 0.5). Defined as
#' \eqn{\rho = \frac{1}{n^2}\sum_{ij} 4 (C_{ij} - 0.5)^2} (Kim & Park, 2007):
#' 1 means clustering is identical across all runs, 0 means maximally unstable.
#'
#' @param consensus consensus matrix from \code{\link{.consensusMatrix}}
#' @return numeric scalar in \code{[0, 1]}
#' @keywords internal
#' @noRd
.dispersionCoefficient <- function(consensus) {
    n <- nrow(consensus)
    sum(4 * (consensus - 0.5)^2) / (n * n)
}


#' Pairwise ARI over all unordered pairs in a list of partitions
#' @param partitions list of integer label vectors
#' @return numeric vector of ARI values, one per pair
#' @keywords internal
#' @noRd
.pairwiseARI <- function(partitions) {
    k <- length(partitions)
    if (k < 2) return(numeric(0))
    pairs <- utils::combn(k, 2)
    vapply(seq_len(ncol(pairs)), function(i) {
        .adjustedRandIndex(partitions[[pairs[1, i]]],
                           partitions[[pairs[2, i]]])
    }, numeric(1))
}


#' Fit one topic model, with or without the PPI/edge view
#'
#' Thin dispatcher used by both \code{compareTopicModels} and
#' \code{bootstrapTopicModels} so the joint-vs-text-only branch lives in one
#' place. Always returns a model carrying \code{W}, \code{H_text}, and
#' \code{H_edges} (folded from the text topics when the edges are excluded).
#'
#' @param X_text paper-by-word matrix
#' @param X_edges paper-by-edge matrix
#' @param k number of topics
#' @param include_ppi logical; jointly factorize the edge view when TRUE
#' @param max_iter,tol,seed passed to the underlying NMF
#' @return fitted model list (W, H_text, H_edges, objective, n_iter)
#' @keywords internal
#' @noRd
.fitTopicModel <- function(X_text, X_edges, k, include_ppi,
                           max_iter = 200, tol = 1e-4, seed = 1) {
    if (include_ppi) {
        model <- .jointNMF(X_text, X_edges, k = k,
                           max_iter = max_iter, tol = tol, seed = seed)
    } else {
        model <- .textNMF(X_text, k = k,
                          max_iter = max_iter, tol = tol, seed = seed)
        model$H_edges <- .edgeLoadingsFromTopics(model$W, X_edges)
    }
    model
}


#' Row-wise cosine similarity between two topic-by-word matrices
#'
#' @param A,B numeric matrices with the same columns (k x p)
#' @return matrix of cosine similarities, \code{nrow(A)} x \code{nrow(B)}
#' @keywords internal
#' @noRd
.cosineSimMatrix <- function(A, B) {
    an <- sqrt(rowSums(A^2)); an[an == 0] <- 1
    bn <- sqrt(rowSums(B^2)); bn[bn == 0] <- 1
    tcrossprod(A, B) / outer(an, bn)
}


#' Greedily match a resample's topics to reference topics
#'
#' Topic indices are arbitrary across NMF fits (label switching), so before any
#' per-topic aggregation each resample's topics must be put back into
#' correspondence with a reference fit. Matches on cosine similarity of the
#' topic-word vectors (\code{H_text} rows), greedily taking the highest
#' remaining similarity until every reference topic has a unique partner. With
#' the small topic counts used here this is effectively optimal.
#'
#' @param sim reference-by-resample similarity matrix from
#'   \code{\link{.cosineSimMatrix}}
#' @return integer vector of length \code{nrow(sim)}: for each reference topic,
#'   the index of the resample topic assigned to it
#' @keywords internal
#' @noRd
.greedyTopicMatch <- function(sim) {
    k <- nrow(sim)
    match <- integer(k)
    used <- logical(ncol(sim))
    remaining <- seq_len(k)
    for (step in seq_len(k)) {
        cand <- which(!used)
        sub <- sim[remaining, cand, drop = FALSE]
        best <- which(sub == max(sub), arr.ind = TRUE)[1, ]
        r <- remaining[best[1]]
        c <- cand[best[2]]
        match[r] <- c
        used[c] <- TRUE
        remaining <- setdiff(remaining, r)
    }
    match
}
