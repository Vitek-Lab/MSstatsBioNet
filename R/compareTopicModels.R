#' Test whether including PPIs changes topic structure beyond random chance
#'
#' Quantifies how much the topic decomposition produced by
#' \code{\link{decomposeSubnetworkByTopic}} changes when the PPI/edge view is
#' included (\code{include_ppi = TRUE}) versus excluded
#' (\code{include_ppi = FALSE}), and separates that change from the run-to-run
#' variability that NMF produces just from its random initialization.
#'
#' NMF converges to a local optimum that depends on the random seed, so a single
#' joint-vs-text comparison conflates the real effect of the PPI view with
#' optimization noise. This function instead refits both modes across many seeds
#' and compares three distributions of partition agreement (Adjusted Rand Index,
#' ARI):
#' \describe{
#'   \item{within_joint}{ARI between pairs of joint runs (different seeds) —
#'     how much the joint solution wobbles on its own.}
#'   \item{within_text}{ARI between pairs of text-only runs — the same for the
#'     text-only solution.}
#'   \item{between}{ARI between the joint and text-only run at the \emph{same}
#'     seed. Because both modes draw \code{W} and \code{H_text} from the same
#'     seeded stream, a matched seed gives both modes an identical
#'     initialization, so this isolates the effect of adding the PPI view from
#'     the starting point.}
#' }
#' If the between-mode ARI is systematically lower than the within-mode ARIs,
#' the PPI view changes the topic structure more than chance would — a
#' one-sided Wilcoxon rank-sum test (\code{between < within}) puts a p-value on
#' it. If the between distribution sits inside the within distributions, the
#' apparent difference is just optimization noise.
#'
#' The expensive, network-bound steps (evidence extraction, abstract fetching,
#' matrix construction) run once; only the NMF is repeated per seed.
#'
#' @param subnetwork list with \code{nodes} and \code{edges} data.frames, e.g.
#'   the output of \code{\link{getSubnetworkFromIndra}}.
#' @param seeds integer vector of NMF seeds to fit (at least 2). Default
#'   \code{1:20}.
#' @param n_topics number of topics (rank of the factorization). Default 5.
#' @param unit either \code{"edges"} (compare edge-to-topic assignments, the
#'   default, matching the subnetworks the decomposition returns) or
#'   \code{"papers"} (compare paper-to-topic assignments).
#' @param min_term_count minimum corpus frequency for a word to be kept when
#'   building the text matrix. Default 2.
#' @param max_iter maximum number of NMF multiplicative-update iterations.
#'   Default 200.
#' @param tol relative-change tolerance for NMF early stopping. Default 1e-4.
#'
#' @return A list with
#'   \describe{
#'     \item{unit}{the comparison unit used.}
#'     \item{seeds}{the seeds fitted.}
#'     \item{n_topics}{the effective number of topics.}
#'     \item{ari}{list of numeric vectors \code{within_joint},
#'       \code{within_text}, and \code{between} (matched seeds).}
#'     \item{summary}{data.frame of median/mean ARI and count per comparison.}
#'     \item{test}{the \code{\link[stats]{wilcox.test}} object comparing the
#'       between distribution against the pooled within distributions
#'       (\code{alternative = "less"}).}
#'     \item{consensus}{list of consensus (co-membership) matrices,
#'       \code{joint} and \code{text}, across seeds.}
#'     \item{dispersion}{named numeric vector of consensus dispersion
#'       coefficients (1 = identical clustering across all seeds).}
#'     \item{partitions}{list of the raw per-seed partitions, \code{joint} and
#'       \code{text}, for further inspection.}
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
#' cmp <- compareTopicModels(subnetwork, seeds = 1:20, n_topics = 5)
#' cmp$summary
#' cmp$test            # p < 0.05 => PPI changes topics beyond chance
#' cmp$dispersion      # how stable each mode is across seeds
#' }
compareTopicModels <- function(subnetwork,
                               seeds = seq_len(20),
                               n_topics = 5,
                               unit = c("edges", "papers"),
                               min_term_count = 2,
                               max_iter = 200,
                               tol = 1e-4) {

    unit <- match.arg(unit)
    # Reuse the decompose validator for subnetwork/n_topics structure.
    .validateDecomposeSubnetworkByTopicInput(subnetwork, n_topics, 0.2)
    if (!is.numeric(seeds) || length(seeds) < 2L || anyNA(seeds) ||
        any(seeds != as.integer(seeds))) {
        stop("`seeds` must be a vector of at least two integer seeds.")
    }
    seeds <- as.integer(seeds)

    # Build the shared matrices once; only the NMF is repeated per seed.
    mats <- .buildTopicMatrices(subnetwork, n_topics, min_term_count)
    X_text <- mats$X_text
    X_edges <- mats$X_edges
    k <- mats$n_topics

    joint_parts <- vector("list", length(seeds))
    text_parts <- vector("list", length(seeds))
    between <- numeric(length(seeds))

    for (i in seq_along(seeds)) {
        s <- seeds[i]
        joint <- .fitTopicModel(X_text, X_edges, k, include_ppi = TRUE,
                                max_iter = max_iter, tol = tol, seed = s)
        text <- .fitTopicModel(X_text, X_edges, k, include_ppi = FALSE,
                               max_iter = max_iter, tol = tol, seed = s)

        joint_parts[[i]] <- .topicPartition(joint, unit)
        text_parts[[i]] <- .topicPartition(text, unit)
        between[i] <- .adjustedRandIndex(joint_parts[[i]], text_parts[[i]])
    }
    names(joint_parts) <- names(text_parts) <- paste0("seed_", seeds)

    within_joint <- .pairwiseARI(joint_parts)
    within_text <- .pairwiseARI(text_parts)
    within_pooled <- c(within_joint, within_text)

    # Is the between-mode agreement lower than within-mode agreement?
    test <- stats::wilcox.test(between, within_pooled,
                               alternative = "less", exact = FALSE)

    consensus_joint <- .consensusMatrix(joint_parts)
    consensus_text <- .consensusMatrix(text_parts)

    summary <- data.frame(
        comparison = c("within_joint", "within_text", "between"),
        n          = c(length(within_joint), length(within_text),
                       length(between)),
        median_ari = c(stats::median(within_joint),
                       stats::median(within_text),
                       stats::median(between)),
        mean_ari   = c(mean(within_joint), mean(within_text), mean(between)),
        stringsAsFactors = FALSE
    )

    list(
        unit       = unit,
        seeds      = seeds,
        n_topics   = k,
        ari        = list(within_joint = within_joint,
                          within_text  = within_text,
                          between       = between),
        summary    = summary,
        test       = test,
        consensus  = list(joint = consensus_joint, text = consensus_text),
        dispersion = c(joint = .dispersionCoefficient(consensus_joint),
                       text  = .dispersionCoefficient(consensus_text)),
        partitions = list(joint = joint_parts, text = text_parts)
    )
}
