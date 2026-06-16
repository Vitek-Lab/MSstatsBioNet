#' Bootstrap the topic decomposition to find each topic's robust top words
#'
#' Refits the NMF topic model on many bootstrap resamples of the papers and
#' reports, for every topic, how reliably each word stays among the topic's top
#' terms. This separates words that genuinely characterise a topic from words
#' that only surface in a single lucky fit, and lets you see how the top-word
#' lists change with and without the PPI view (run once per \code{include_ppi}).
#'
#' Because topic indices are arbitrary across fits (label switching), each
#' resample's topics are first aligned to a reference fit on the full data by
#' cosine similarity of their topic-word vectors. The papers are resampled with
#' replacement; the word vocabulary is held fixed (from the full data) so topics
#' remain comparable across resamples, and the NMF seed is held fixed so the
#' variability reported reflects \emph{data} resampling rather than random
#' initialization.
#'
#' @param subnetwork list with \code{nodes} and \code{edges} data.frames, e.g.
#'   the output of \code{\link{getSubnetworkFromIndra}}.
#' @param n_boot number of bootstrap resamples. Default 50.
#' @param n_topics number of topics (rank of the factorization). Default 5.
#' @param include_ppi logical; factorize the PPI/edge view jointly with the text
#'   (\code{TRUE}, default) or use paper words only (\code{FALSE}). See
#'   \code{\link{decomposeSubnetworkByTopic}}.
#' @param n_top_terms number of top words that define a topic's "top list" in
#'   each resample (the cutoff for the selection-frequency tally). Default 10.
#' @param min_term_count minimum corpus frequency for a word to be kept when
#'   building the text matrix. Default 2.
#' @param max_iter maximum number of NMF multiplicative-update iterations.
#'   Default 200.
#' @param tol relative-change tolerance for NMF early stopping. Default 1e-4.
#' @param seed random seed for the reference fit, the resampling, and each
#'   bootstrap NMF. Default 1.
#'
#' @return A list with
#'   \describe{
#'     \item{include_ppi, n_boot, n_topics}{the settings used.}
#'     \item{topTerms}{named list \code{topic_1} ... \code{topic_k}. Each is a
#'       data.frame sorted by \code{selection_freq}, with columns \code{term},
#'       \code{selection_freq} (fraction of resamples the word was in this
#'       topic's top \code{n_top_terms}), and \code{mean_weight} (mean
#'       within-topic word weight across resamples). A word with
#'       \code{selection_freq} near 1 is a stable signature of the topic.}
#'     \item{reference}{named list of the top \code{n_top_terms} words per topic
#'       from the single full-data fit, for comparison.}
#'   }
#'
#' @seealso \code{\link{decomposeSubnetworkByTopic}},
#'   \code{\link{compareTopicModels}}
#' @note **Beta feature:** This function is experimental and the API may
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
#'
#' # Top words with PPIs included vs. words only:
#' boot_ppi  <- bootstrapTopicModels(subnetwork, include_ppi = TRUE)
#' boot_text <- bootstrapTopicModels(subnetwork, include_ppi = FALSE)
#'
#' head(boot_ppi$topTerms$topic_1)     # robust signature words for topic 1
#' boot_text$topTerms$topic_1
#' }
bootstrapTopicModels <- function(subnetwork,
                                 n_boot = 50,
                                 n_topics = 5,
                                 include_ppi = TRUE,
                                 n_top_terms = 10,
                                 min_term_count = 2,
                                 max_iter = 200,
                                 tol = 1e-4,
                                 seed = 1) {

    .validateDecomposeSubnetworkByTopicInput(subnetwork, n_topics, 0.2,
                                             include_ppi)
    if (!is.numeric(n_boot) || length(n_boot) != 1L || is.na(n_boot) ||
        n_boot < 2 || n_boot != as.integer(n_boot)) {
        stop("`n_boot` must be a single integer >= 2.")
    }
    n_boot <- as.integer(n_boot)

    # Build the shared matrices once; keep the vocabulary fixed across resamples.
    mats <- .buildTopicMatrices(subnetwork, n_topics, min_term_count)
    X_text <- mats$X_text
    X_edges <- mats$X_edges
    k <- mats$n_topics
    n_papers <- nrow(X_text)
    terms <- colnames(X_text)

    # Reference fit on the full data defines the canonical topic ordering that
    # every resample is aligned back to.
    ref <- .fitTopicModel(X_text, X_edges, k, include_ppi,
                          max_iter = max_iter, tol = tol, seed = seed)
    reference <- lapply(seq_len(k), function(t)
        .topTermsForTopic(ref$H_text, t, n_top_terms))
    names(reference) <- paste0("topic_", seq_len(k))

    # Precompute all resample row-index sets up front: the NMF resets the RNG
    # internally, so drawing the samples inside the loop would repeat them.
    set.seed(seed)
    boot_indices <- lapply(seq_len(n_boot),
                           function(b) sample.int(n_papers, replace = TRUE))

    # Tally, per (reference topic x word): top-list selection count and the
    # summed within-topic weight, accumulated across resamples.
    topcount <- matrix(0, k, length(terms), dimnames = list(NULL, terms))
    weightsum <- matrix(0, k, length(terms), dimnames = list(NULL, terms))

    for (b in seq_len(n_boot)) {
        idx <- boot_indices[[b]]
        fit <- .fitTopicModel(X_text[idx, , drop = FALSE],
                              X_edges[idx, , drop = FALSE],
                              k, include_ppi,
                              max_iter = max_iter, tol = tol, seed = seed)

        # Align this resample's topics to the reference, then aggregate.
        map <- .greedyTopicMatch(.cosineSimMatrix(ref$H_text, fit$H_text))
        for (r in seq_len(k)) {
            row <- fit$H_text[map[r], ]
            s <- sum(row)
            weightsum[r, ] <- weightsum[r, ] + (if (s > 0) row / s else row)
            top <- order(row, decreasing = TRUE)[seq_len(n_top_terms)]
            top <- top[row[top] > 0]
            topcount[r, top] <- topcount[r, top] + 1
        }
    }

    topTerms <- lapply(seq_len(k), function(r) {
        keep <- which(topcount[r, ] > 0)
        df <- data.frame(
            term            = terms[keep],
            selection_freq  = topcount[r, keep] / n_boot,
            mean_weight     = weightsum[r, keep] / n_boot,
            stringsAsFactors = FALSE
        )
        df <- df[order(-df$selection_freq, -df$mean_weight), , drop = FALSE]
        rownames(df) <- NULL
        df
    })
    names(topTerms) <- paste0("topic_", seq_len(k))

    list(
        include_ppi = include_ppi,
        n_boot      = n_boot,
        n_topics    = k,
        topTerms    = topTerms,
        reference   = reference
    )
}
