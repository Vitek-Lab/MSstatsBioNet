# Synthetic corpus: three themes, each with its own vocabulary, papers, and
# edges, so the NMF has clear structure to recover without any network calls.
THEME_WORDS <- list(
    c("kinase", "phosphorylation", "signaling", "cascade", "mapk", "erk"),
    c("dna", "repair", "damage", "checkpoint", "replication", "genome"),
    c("immune", "cytokine", "inflammation", "macrophage", "tcell", "interferon")
)

make_theme_corpus <- function(edges_per_theme = 15, papers_per_theme = 8) {
    edges <- list()
    evidence <- list()
    abstracts <- character(0)
    for (th in seq_along(THEME_WORDS)) {
        pmids <- paste0("PM", th, "_", seq_len(papers_per_theme))
        words <- THEME_WORDS[[th]]
        abstracts[pmids] <- vapply(seq_along(pmids), function(i) {
            paste(rep(words, times = 3 + (i %% 3)), collapse = " ")
        }, character(1))
        for (e in seq_len(edges_per_theme)) {
            src <- paste0("G", th, "_", e)
            tgt <- paste0("G", th, "_", e + 1)
            hash <- paste0("h", th, "_", e)
            edges[[length(edges) + 1]] <- data.frame(
                source = src, target = tgt, interaction = "Activation",
                site = NA_character_, evidenceLink = "https://example.com",
                stmt_hash = hash, stringsAsFactors = FALSE
            )
            ev_pmids <- pmids[c(e %% papers_per_theme + 1,
                                (e + 3) %% papers_per_theme + 1)]
            evidence[[length(evidence) + 1]] <- data.frame(
                source = src, target = tgt, interaction = "Activation",
                site = NA_character_, evidenceLink = "https://example.com",
                stmt_hash = hash, text = "sentence", pmid = ev_pmids,
                stringsAsFactors = FALSE
            )
        }
    }
    edges <- do.call(rbind, edges)
    nodes <- data.frame(id = unique(c(edges$source, edges$target)),
                        stringsAsFactors = FALSE)
    list(subnetwork = list(nodes = nodes, edges = edges),
         evidence = do.call(rbind, evidence),
         abstracts = abstracts)
}

forbid_network <- function(env = parent.frame()) {
    testthat::local_mocked_bindings(
        .extract_evidence_text = function(...) stop("INDRA was queried"),
        .fetch_clean_abstracts_xml = function(...) stop("PubMed was queried"),
        .env = env
    )
}

describe("decomposeSubnetworkByTopic with a supplied corpus", {

    test_that("uses the supplied evidence and abstracts without fetching", {
        forbid_network()
        corpus <- make_theme_corpus()
        topics <- decomposeSubnetworkByTopic(
            corpus$subnetwork, n_topics = 3,
            evidence = corpus$evidence, abstracts = corpus$abstracts
        )
        expect_length(topics, 3)
        used <- attr(topics, "corpus")
        expect_setequal(names(used$abstracts), unique(corpus$evidence$pmid))
        expect_equal(nrow(used$evidence), nrow(corpus$evidence))
    })

    test_that("a topic subnetwork can be re-decomposed from the corpus attr", {
        forbid_network()
        corpus <- make_theme_corpus()
        topics <- decomposeSubnetworkByTopic(
            corpus$subnetwork, n_topics = 3, edge_topic_cutoff = 0.9,
            evidence = corpus$evidence, abstracts = corpus$abstracts
        )
        used <- attr(topics, "corpus")
        deeper <- decomposeSubnetworkByTopic(
            topics$topic_1, n_topics = 2,
            evidence = used$evidence, abstracts = used$abstracts
        )
        # Evidence is restricted to the parent topic's edges.
        sub_used <- attr(deeper, "corpus")$evidence
        expect_true(all(sub_used$stmt_hash %in% topics$topic_1$edges$stmt_hash))
    })

    test_that("only PMIDs missing from `abstracts` are fetched", {
        corpus <- make_theme_corpus()
        missing <- names(corpus$abstracts)[1:2]
        fetched <- NULL
        testthat::local_mocked_bindings(
            .fetch_clean_abstracts_xml = function(pmids, ...) {
                fetched <<- pmids
                as.list(stats::setNames(corpus$abstracts[pmids], pmids))
            }
        )
        decomposeSubnetworkByTopic(
            corpus$subnetwork, n_topics = 3, evidence = corpus$evidence,
            abstracts = corpus$abstracts[-(1:2)]
        )
        expect_setequal(fetched, missing)
    })

    test_that("rejects malformed evidence and abstracts", {
        corpus <- make_theme_corpus()
        expect_error(
            decomposeSubnetworkByTopic(corpus$subnetwork,
                                       evidence = data.frame(x = 1)),
            "`evidence` must be a data.frame"
        )
        expect_error(
            decomposeSubnetworkByTopic(corpus$subnetwork,
                                       evidence = corpus$evidence,
                                       abstracts = unname(corpus$abstracts)),
            "`abstracts` must be a named"
        )
    })
})

describe("decomposeSubnetworkIntoHierarchicalTopics", {

    test_that("splits until every leaf is small or cannot split further", {
        forbid_network()
        corpus <- make_theme_corpus()
        h <- decomposeSubnetworkIntoHierarchicalTopics(
            corpus$subnetwork, max_edges = 10, n_topics = 3,
            evidence = corpus$evidence, abstracts = corpus$abstracts
        )
        expect_s3_class(h, "topicHierarchy")
        tree <- h$tree
        expect_equal(tree$id[1], "root")
        expect_true(is.na(tree$parent_id[1]))
        expect_equal(tree$n_edges[1], nrow(corpus$subnetwork$edges))
        expect_true(any(tree$depth > 0))

        leaves <- tree[tree$is_leaf, ]
        expect_true(all(leaves$n_edges <= 10 |
                            leaves$stop_reason != "small_enough"))
        expect_true(all(is.na(tree$stop_reason[!tree$is_leaf])))
    })

    test_that("tree, subnetworks, and edge membership are consistent", {
        forbid_network()
        corpus <- make_theme_corpus()
        h <- decomposeSubnetworkIntoHierarchicalTopics(
            corpus$subnetwork, max_edges = 10, n_topics = 3,
            evidence = corpus$evidence, abstracts = corpus$abstracts
        )
        tree <- h$tree
        expect_setequal(names(h$subnetworks), tree$id)
        expect_true(all(tree$parent_id[-1] %in% tree$id))
        expect_equal(unname(vapply(h$subnetworks[tree$id],
                                   function(s) nrow(s$edges), integer(1))),
                     tree$n_edges)

        # Children are strict subsets of their parent's edges.
        for (i in which(!is.na(tree$parent_id))) {
            child <- h$subnetworks[[tree$id[i]]]$edges
            parent <- h$subnetworks[[tree$parent_id[i]]]$edges
            expect_true(all(child$stmt_hash %in% parent$stmt_hash))
            expect_lt(nrow(child), nrow(parent))
        }

        expect_equal(nrow(h$edge_membership), sum(tree$n_edges[-1]))
        expect_equal(h$tree$pathString[1], "root")
        child_row <- which(tree$depth == 1)[1]
        expect_equal(tree$pathString[child_row],
                     paste0("root/", tree$id[child_row]))
    })

    test_that("queries INDRA and PubMed once when no corpus is supplied", {
        corpus <- make_theme_corpus()
        calls <- c(indra = 0, pubmed = 0)
        testthat::local_mocked_bindings(
            .extract_evidence_text = function(df) {
                calls[["indra"]] <<- calls[["indra"]] + 1
                corpus$evidence
            },
            .fetch_clean_abstracts_xml = function(pmids, ...) {
                calls[["pubmed"]] <<- calls[["pubmed"]] + 1
                as.list(stats::setNames(corpus$abstracts[pmids], pmids))
            }
        )
        h <- decomposeSubnetworkIntoHierarchicalTopics(
            corpus$subnetwork, max_edges = 10, n_topics = 3
        )
        expect_gt(nrow(h$tree), 1)
        expect_equal(calls, c(indra = 1, pubmed = 1))
    })

    test_that("max_depth = 0 returns only the root", {
        forbid_network()
        corpus <- make_theme_corpus()
        h <- decomposeSubnetworkIntoHierarchicalTopics(
            corpus$subnetwork, max_depth = 0,
            evidence = corpus$evidence, abstracts = corpus$abstracts
        )
        expect_equal(nrow(h$tree), 1)
        expect_equal(h$tree$stop_reason, "max_depth")
        expect_equal(nrow(h$edge_membership), 0)
    })

    test_that("a small network is a single small_enough leaf", {
        forbid_network()
        corpus <- make_theme_corpus(edges_per_theme = 3)
        h <- decomposeSubnetworkIntoHierarchicalTopics(
            corpus$subnetwork, max_edges = 10,
            evidence = corpus$evidence, abstracts = corpus$abstracts
        )
        expect_equal(nrow(h$tree), 1)
        expect_equal(h$tree$stop_reason, "small_enough")
    })

    test_that("print shows an indented tree", {
        forbid_network()
        corpus <- make_theme_corpus()
        h <- decomposeSubnetworkIntoHierarchicalTopics(
            corpus$subnetwork, max_edges = 10, n_topics = 3,
            evidence = corpus$evidence, abstracts = corpus$abstracts
        )
        out <- capture.output(print(h))
        expect_match(out[1], "^Topic hierarchy:")
        expect_match(out[2], "^root \\[45 edges")
        expect_true(any(grepl("^  1 \\[", out)))
    })

    test_that("validates max_edges and max_depth", {
        corpus <- make_theme_corpus()
        expect_error(
            decomposeSubnetworkIntoHierarchicalTopics(corpus$subnetwork,
                                                      max_edges = 0),
            "`max_edges`"
        )
        expect_error(
            decomposeSubnetworkIntoHierarchicalTopics(corpus$subnetwork,
                                                      max_depth = 1.5),
            "`max_depth`"
        )
    })
})
