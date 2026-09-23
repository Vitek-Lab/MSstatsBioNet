make_edges <- function() {
    data.frame(
        source       = c("A", "B"),
        target       = c("B", "C"),
        interaction  = c("activates", "inhibits"),
        site         = c("T308", "S473"),
        evidenceLink = c("https://example.com/1", "https://example.com/2"),
        stmt_hash    = c("hash1", "hash2"),
        stringsAsFactors = FALSE
    )
}

make_nodes <- function() {
    data.frame(
        id    = c("A", "B", "C"),
        label = c("GeneA", "GeneB", "GeneC"),
        stringsAsFactors = FALSE
    )
}

make_mock_evidence <- function() {
    data.frame(
        source       = c("A", "B"),
        target       = c("B", "C"),
        interaction  = c("activates", "inhibits"),
        site         = c("T308", "S473"),
        evidenceLink = c("https://example.com/1", "https://example.com/2"),
        stmt_hash    = c("hash1", "hash2"),
        text         = c("CHEK1 sentence.", "Lipid sentence."),
        pmid         = c("11111111", "22222222"),
        stringsAsFactors = FALSE
    )
}

make_mock_abstracts <- function() {
    list(
        "11111111" = "CHEK1 phosphorylates CDC25A in response to DNA damage.",
        "22222222" = "Unrelated text about lipid metabolism and glucose uptake."
    )
}

make_pubmed_xml <- function(pmids) {
    articles <- vapply(pmids, function(pmid) {
        sprintf(paste0(
            "<PubmedArticle><MedlineCitation><PMID>%s</PMID><Article><Abstract>",
            "<AbstractText>Abstract for %s.</AbstractText>",
            "</Abstract></Article></MedlineCitation></PubmedArticle>"
        ), pmid, pmid)
    }, character(1))

    paste0("<PubmedArticleSet>", paste(articles, collapse = ""), "</PubmedArticleSet>")
}

describe(".score_by_tag_count", {
    
    test_that("returns 0 for an abstract that contains none of the tags", {
        scores <- .score_by_tag_count("nothing relevant here", c("CHEK1", "DNA damage"))
        expect_equal(scores, 0L)
    })
    
    test_that("counts multiple matching tags correctly", {
        abstract <- "CHEK1 is involved in DNA damage repair pathways."
        scores   <- .score_by_tag_count(abstract, c("chek1", "DNA damage", "apoptosis"))
        expect_equal(scores, 2L)
    })
    
})

describe(".contains_any_keyword", {

    test_that("matches whole words only, not longer words that contain the keyword", {
        abstracts <- c(
            "Tumors of the colon were resected.",
            "Bacteria colonize the gut; colonization was measured.",
            "Colon-specific expression was high.",
            "COLON cancer cohort."
        )
        expect_equal(.contains_any_keyword(abstracts, "colon"),
                     c(TRUE, FALSE, TRUE, TRUE))
    })

    test_that("matches multi-word phrases and treats regex characters literally", {
        abstracts <- c("Levels of IL-6 (pg/ml) rose.", "DNA damage repair.",
                       "ILx6 was absent.")
        expect_equal(.contains_any_keyword(abstracts, c("il-6", "dna damage")),
                     c(TRUE, TRUE, FALSE))
        expect_equal(.contains_any_keyword(abstracts, "(pg/ml)"),
                     c(TRUE, FALSE, FALSE))
    })

})

describe(".score_by_cosine", {
    
    test_that("returns a numeric vector of the same length as abstracts", {
        abstracts <- c("DNA damage repair involves CHEK1.", "Unrelated text about metabolism.")
        scores    <- .score_by_cosine("CHEK1 DNA damage", abstracts)
        expect_true(is.numeric(scores))
        expect_length(scores, 2)
    })
    
    test_that("scores a highly relevant abstract higher than an irrelevant one", {
        relevant   <- "CHEK1 mediates the DNA damage checkpoint response."
        irrelevant <- "Photosynthesis occurs in the chloroplast of plant cells."
        scores     <- .score_by_cosine("CHEK1 DNA damage checkpoint", c(relevant, irrelevant))
        expect_gt(scores[1], scores[2])
    })
    
})

describe(".extract_evidence_text", {
    
    test_that("stops when required columns are missing from the edges dataframe", {
        bad_df <- data.frame(source = "A", target = "B", stringsAsFactors = FALSE)
        expect_error(
            .extract_evidence_text(bad_df),
            regexp = "Missing required columns"
        )
    })
    
    test_that("returns an empty dataframe with correct columns when the INDRA API returns nothing", {
        edges <- make_edges()
        
        # Mock .query_indra_evidence to always return NULL
        mockery::stub(.extract_evidence_text, ".query_indra_evidence", NULL)
        
        result <- suppressWarnings(.extract_evidence_text(edges))
        expect_s3_class(result, "data.frame")
        expect_true(all(c("source", "target", "interaction", "site",
                          "evidenceLink", "stmt_hash", "text", "pmid") %in% names(result)))
        expect_equal(nrow(result), 0)
    })
    
})

describe(".fetch_clean_abstracts_xml", {
    
    test_that("returns an empty list when given an empty pmids vector", {
        result <- .fetch_clean_abstracts_xml(character(0))
        expect_true(is.list(result))
        expect_length(result, 0)
    })
    
    test_that("stores an empty string for a PMID that triggers an API error", {
        mockery::stub(
            .fetch_clean_abstracts_xml,
            "entrez_fetch",
            function(...) stop("network error")
        )
        result <- suppressMessages(.fetch_clean_abstracts_xml(c("99999999")))
        expect_true("99999999" %in% names(result))
        expect_equal(result[["99999999"]], "")
    })

    test_that("requests PMIDs in batches and keeps every PMID in the result", {
        requested <- list()
        mockery::stub(
            .fetch_clean_abstracts_xml,
            "entrez_fetch",
            function(db, id, rettype) {
                requested[[length(requested) + 1]] <<- id
                make_pubmed_xml(id)
            }
        )
        pmids  <- c("11111111", "11111112", "11111113", "11111114", "11111115")
        result <- suppressMessages(
            .fetch_clean_abstracts_xml(pmids, batch_size = 2)
        )

        expect_equal(
            requested,
            list(c("11111111", "11111112"), c("11111113", "11111114"), "11111115")
        )
        expect_equal(names(result), pmids)
        expect_equal(result[["11111113"]], "Abstract for 11111113.")
    })

    test_that("returns an empty string for a PMID missing from the response", {
        mockery::stub(
            .fetch_clean_abstracts_xml,
            "entrez_fetch",
            function(db, id, rettype) make_pubmed_xml(setdiff(id, "22222222"))
        )
        result <- suppressMessages(
            .fetch_clean_abstracts_xml(c("11111111", "22222222"))
        )
        expect_equal(result[["11111111"]], "Abstract for 11111111.")
        expect_equal(result[["22222222"]], "")
    })

})

describe(".parse_pubmed_abstracts", {

    test_that("joins multi-section abstracts and ignores cited reference PMIDs", {
        record <- paste0(
            "<PubmedArticleSet><PubmedArticle>",
            "<MedlineCitation><PMID>11111111</PMID>",
            "<Article><Abstract>",
            "<AbstractText Label='BACKGROUND'>  First part. </AbstractText>",
            "<AbstractText Label='RESULTS'>Second part.</AbstractText>",
            "</Abstract></Article></MedlineCitation>",
            "<PubmedData><ReferenceList><Reference><ArticleIdList>",
            "<ArticleId IdType='pubmed'>99999999</ArticleId>",
            "</ArticleIdList></Reference></ReferenceList></PubmedData>",
            "</PubmedArticle></PubmedArticleSet>"
        )
        result <- .parse_pubmed_abstracts(record)
        expect_equal(names(result), "11111111")
        expect_equal(result[["11111111"]], "First part. Second part.")
    })

    test_that("returns an empty string for an article with no abstract", {
        record <- paste0(
            "<PubmedArticleSet><PubmedArticle><MedlineCitation>",
            "<PMID>11111111</PMID><Article><ArticleTitle>T</ArticleTitle></Article>",
            "</MedlineCitation></PubmedArticle></PubmedArticleSet>"
        )
        expect_equal(.parse_pubmed_abstracts(record)[["11111111"]], "")
    })

    test_that("returns an empty list when the response has no articles", {
        result <- .parse_pubmed_abstracts("<PubmedArticleSet></PubmedArticleSet>")
        expect_true(is.list(result))
        expect_length(result, 0)
    })

})

describe("filterSubnetworkByContext", {
    
    test_that("returns nodes, edges, and evidence that match the query tags (happy path)", {
        nodes <- make_nodes()
        edges <- make_edges()
        
        # --- mock .extract_evidence_text ---
        mock_evidence <- data.frame(
            source       = c("A", "B"),
            target       = c("B", "C"),
            interaction  = c("activates", "inhibits"),
            site         = c("T308", "S473"),
            evidenceLink = c("https://example.com/1", "https://example.com/2"),
            stmt_hash    = c("hash1", "hash2"),
            text         = c(
                "CHEK1 phosphorylates CDC25A in response to DNA damage.",
                "Unrelated text about lipid metabolism and glucose uptake."
            ),
            pmid         = c("11111111", "22222222"),
            stringsAsFactors = FALSE
        )
        mockery::stub(filterSubnetworkByContext, ".extract_evidence_text", mock_evidence)
        
        # --- mock .fetch_clean_abstracts_xml ---
        mock_abstracts <- list(
            "11111111" = "CHEK1 phosphorylates CDC25A in response to DNA damage.",
            "22222222" = "Unrelated text about lipid metabolism and glucose uptake."
        )
        mockery::stub(filterSubnetworkByContext, ".fetch_clean_abstracts_xml", mock_abstracts)
        
        result <- filterSubnetworkByContext(
            nodes  = nodes,
            edges  = edges,
            query  = c("CHEK1", "DNA damage"),
            cutoff = 1,
            method = "tag_count"
        )
        
        # Structure check
        expect_named(result, c("nodes", "edges", "evidence", "abstracts"))
        
        # Only the CHEK1/DNA-damage abstract passed the cutoff
        expect_equal(nrow(result$edges), 1)
        expect_equal(result$edges$stmt_hash, "hash1")
        
        # Nodes should only contain those referenced by surviving edges
        expect_true(all(result$nodes$id %in% c("A", "B", "C")))
        expect_false("C" %in% result$nodes$id)   # C only appears in the filtered-out edge
        
        # Evidence rows carry a score column
        expect_true("score" %in% names(result$evidence))
        expect_true(all(result$evidence$score >= 1))
    })

    test_that("returns the abstracts of the PMIDs in the kept evidence", {
        mockery::stub(filterSubnetworkByContext, ".extract_evidence_text",
                      make_mock_evidence())
        mockery::stub(filterSubnetworkByContext, ".fetch_clean_abstracts_xml",
                      make_mock_abstracts())

        result <- filterSubnetworkByContext(
            make_nodes(), make_edges(), query = c("CHEK1", "DNA damage")
        )

        expect_type(result$abstracts, "character")
        expect_equal(names(result$abstracts), "11111111")
        expect_equal(result$abstracts[["11111111"]],
                     make_mock_abstracts()[["11111111"]])
    })

    test_that("exclude_keywords drops abstracts containing any keyword", {
        mockery::stub(filterSubnetworkByContext, ".extract_evidence_text",
                      make_mock_evidence())
        mockery::stub(filterSubnetworkByContext, ".fetch_clean_abstracts_xml",
                      make_mock_abstracts())

        # cutoff = 0 keeps every abstract, so only the exclusion filters.
        result <- filterSubnetworkByContext(
            make_nodes(), make_edges(), query = "CHEK1", cutoff = 0,
            exclude_keywords = c("apoptosis", "LIPID")
        )

        expect_equal(result$edges$stmt_hash, "hash1")
        expect_equal(unique(result$evidence$pmid), "11111111")
        expect_equal(names(result$abstracts), "11111111")
        expect_false("C" %in% result$nodes$id)
    })

    test_that("exclude_keywords takes precedence over a passing score", {
        mockery::stub(filterSubnetworkByContext, ".extract_evidence_text",
                      make_mock_evidence())
        mockery::stub(filterSubnetworkByContext, ".fetch_clean_abstracts_xml",
                      make_mock_abstracts())

        result <- filterSubnetworkByContext(
            make_nodes(), make_edges(), query = "CHEK1", cutoff = 1,
            exclude_keywords = "cdc25a"
        )

        expect_equal(nrow(result$edges), 0)
        expect_equal(nrow(result$evidence), 0)
        expect_length(result$abstracts, 0)
    })

    test_that("exclude_keywords works with the cosine method", {
        mockery::stub(filterSubnetworkByContext, ".extract_evidence_text",
                      make_mock_evidence())
        mockery::stub(filterSubnetworkByContext, ".fetch_clean_abstracts_xml",
                      make_mock_abstracts())

        result <- filterSubnetworkByContext(
            make_nodes(), make_edges(), query = "CHEK1 DNA damage",
            cutoff = 0, method = "cosine", exclude_keywords = "glucose"
        )

        expect_equal(result$edges$stmt_hash, "hash1")
    })

    test_that("filters by exclude_keywords alone when query is omitted", {
        mockery::stub(filterSubnetworkByContext, ".extract_evidence_text",
                      make_mock_evidence())
        mockery::stub(filterSubnetworkByContext, ".fetch_clean_abstracts_xml",
                      make_mock_abstracts())

        result <- filterSubnetworkByContext(
            make_nodes(), make_edges(), exclude_keywords = "lipid"
        )

        expect_equal(result$edges$stmt_hash, "hash1")
        expect_equal(names(result$abstracts), "11111111")
        expect_true(all(is.na(result$evidence$score)))
    })

    test_that("keeps every abstract when no exclude keyword matches", {
        mockery::stub(filterSubnetworkByContext, ".extract_evidence_text",
                      make_mock_evidence())
        mockery::stub(filterSubnetworkByContext, ".fetch_clean_abstracts_xml",
                      make_mock_abstracts())

        result <- filterSubnetworkByContext(
            make_nodes(), make_edges(), exclude_keywords = "photosynthesis"
        )

        expect_equal(nrow(result$edges), 2)
        expect_setequal(names(result$abstracts), c("11111111", "22222222"))
    })

    test_that("requires query or exclude_keywords, and cutoff needs query", {
        expect_error(
            filterSubnetworkByContext(make_nodes(), make_edges()),
            "Supply `query`, `exclude_keywords`, or both"
        )
        expect_error(
            filterSubnetworkByContext(make_nodes(), make_edges(), cutoff = 1,
                                      exclude_keywords = "lipid"),
            "`cutoff` has no effect without `query`"
        )
    })

    test_that("rejects malformed exclude_keywords", {
        for (bad in list(1, NA_character_, "", character(0))) {
            expect_error(
                filterSubnetworkByContext(make_nodes(), make_edges(),
                                          query = "CHEK1",
                                          exclude_keywords = bad),
                "`exclude_keywords` must be NULL"
            )
        }
    })

    test_that("returns empty abstracts when no evidence is found", {
        mockery::stub(filterSubnetworkByContext, ".extract_evidence_text",
                      make_mock_evidence()[0, ])

        result <- suppressWarnings(filterSubnetworkByContext(
            make_nodes(), make_edges(), query = "CHEK1"
        ))

        expect_named(result, c("nodes", "edges", "evidence", "abstracts"))
        expect_length(result$abstracts, 0)
    })

    test_that("output can be decomposed without re-querying INDRA or PubMed", {
        themes <- list(
            c("kinase", "phosphorylation", "signaling", "cascade", "mapk"),
            c("dna", "repair", "damage", "checkpoint", "replication"),
            c("immune", "cytokine", "inflammation", "macrophage", "interferon")
        )
        edges <- list()
        evidence <- list()
        abstracts <- list()
        for (th in seq_along(themes)) {
            pmids <- paste0(th, "000", seq_len(6))
            for (i in seq_along(pmids)) {
                abstracts[[pmids[i]]] <- paste(
                    rep(themes[[th]], times = 2 + i %% 3), collapse = " "
                )
            }
            for (e in seq_len(8)) {
                edge <- data.frame(
                    source = paste0("G", th, "_", e),
                    target = paste0("G", th, "_", e + 1),
                    interaction = "Activation", site = NA_character_,
                    evidenceLink = "https://example.com",
                    stmt_hash = paste0("h", th, "_", e),
                    stringsAsFactors = FALSE
                )
                edges[[length(edges) + 1]] <- edge
                evidence[[length(evidence) + 1]] <- cbind(
                    edge, text = "sentence",
                    pmid = pmids[c(e %% 6 + 1, (e + 2) %% 6 + 1)],
                    stringsAsFactors = FALSE
                )
            }
        }
        edges <- do.call(rbind, edges)
        nodes <- data.frame(id = unique(c(edges$source, edges$target)),
                            stringsAsFactors = FALSE)
        mockery::stub(filterSubnetworkByContext, ".extract_evidence_text",
                      do.call(rbind, evidence))
        mockery::stub(filterSubnetworkByContext, ".fetch_clean_abstracts_xml",
                      abstracts)

        filtered <- filterSubnetworkByContext(
            nodes, edges, exclude_keywords = "cytokine"
        )
        expect_false(any(grepl("^G3_", filtered$edges$source)))

        testthat::local_mocked_bindings(
            .extract_evidence_text = function(...) stop("INDRA was queried"),
            .fetch_clean_abstracts_xml = function(...) stop("PubMed was queried")
        )
        hierarchy <- decomposeSubnetworkIntoHierarchicalTopics(
            filtered, max_edges = 8, n_topics = 2,
            evidence = filtered$evidence, abstracts = filtered$abstracts
        )

        expect_s3_class(hierarchy, "topicHierarchy")
        expect_equal(hierarchy$tree$n_edges[1], nrow(filtered$edges))
        expect_setequal(names(hierarchy$corpus$abstracts),
                        unique(filtered$evidence$pmid))
        expect_true(any(!hierarchy$tree$is_leaf))
    })

})
