#' Filter a subnetwork by contextual relevance using TF-IDF cosine similarity
#'
#' Fetches PubMed abstracts for evidence PMIDs, scores each abstract against a
#' user-supplied text query, and returns only the nodes, edges, and evidence
#' rows whose abstracts meet the similarity cutoff.
#'
#' @param nodes     A dataframe of network nodes.
#' @param edges     A dataframe of network edges with columns: source, target,
#'                  interaction, site, evidenceLink, stmt_hash.
#' @param similarity_cutoff Numeric in [-1, 1]. Only evidence whose abstract
#'                  scores >= this value is retained. Default 0.10.
#' @param query     Character string. The text query to compare against abstracts.
#'                  Expand with synonyms / related terms for better recall.
#'
#' @return A named list with three elements:
#'   \item{nodes}{Filtered nodes dataframe (only nodes present in kept edges)}
#'   \item{edges}{Filtered edges dataframe}
#'   \item{evidence}{Dataframe with columns: source, target, interaction, site,
#'         evidenceLink, stmt_hash, text, pmid, similarity}
#'
#' @importFrom text2vec itoken word_tokenizer create_vocabulary prune_vocabulary
#'   vocab_vectorizer create_dtm TfIdf fit_transform
#' @importFrom stopwords stopwords
#' @export
filterSubnetworkByContext <- function(nodes,
                                      edges,
                                      similarity_cutoff = 0.10,
                                      query = "DNA damage repair cancer oncology") {
    
    # ── 1. Extract evidence text from edges ───────────────────────────────────
    evidence <- .extract_evidence_text(edges)
    
    if (nrow(evidence) == 0) {
        warning("No evidence text found — returning unfiltered inputs.")
        return(list(nodes = nodes, edges = edges, evidence = evidence))
    }
    
    # ── 2. Fetch PubMed abstracts for unique PMIDs ────────────────────────────
    pmids <- unique(evidence$pmid[nchar(evidence$pmid) > 0])
    
    if (length(pmids) == 0) {
        warning("No PMIDs found in evidence — returning unfiltered inputs.")
        return(list(nodes = nodes, edges = edges, evidence = evidence))
    }
    
    abstract_list   <- .fetch_clean_abstracts_xml(pmids)
    abstracts_df    <- data.frame(
        pmid     = names(abstract_list),
        abstract = unlist(abstract_list, use.names = FALSE),
        stringsAsFactors = FALSE
    )
    
    # ── 3. TF-IDF vectorisation (query + all abstracts) ───────────────────────
    all_texts <- c(query, abstracts_df$abstract)
    
    tokens     <- itoken(all_texts,
                         preprocessor = tolower,
                         tokenizer    = word_tokenizer)
    vocab      <- create_vocabulary(tokens, stopwords = stopwords("en"))
    vocab      <- prune_vocabulary(vocab, term_count_min = 1)
    vectorizer <- vocab_vectorizer(vocab)
    dtm        <- create_dtm(tokens, vectorizer)
    tfidf      <- TfIdf$new()
    dtm_tfidf  <- fit_transform(dtm, tfidf)
    
    # ── 4. Cosine similarity: query (row 1) vs each abstract ─────────────────
    .cos_sim <- function(a, b) {
        a <- as.numeric(a)
        b <- as.numeric(b)
        sum(a * b) / (sqrt(sum(a^2)) * sqrt(sum(b^2)))
    }
    
    query_vec     <- dtm_tfidf[1, , drop = FALSE]
    abstract_vecs <- dtm_tfidf[-1, , drop = FALSE]
    
    scores <- sapply(seq_len(nrow(abstract_vecs)), function(i) {
        .cos_sim(query_vec, abstract_vecs[i, , drop = FALSE])
    })
    
    abstracts_df$similarity <- round(scores, 4)
    
    # ── 5. Filter abstracts by similarity cutoff ──────────────────────────────
    passing_pmids <- abstracts_df$pmid[abstracts_df$similarity >= similarity_cutoff]
    
    cat(sprintf(
        "\n%d / %d abstracts passed similarity cutoff (>= %.2f)\n",
        length(passing_pmids), nrow(abstracts_df), similarity_cutoff
    ))
    
    # ── 6. Filter evidence, edges, nodes ─────────────────────────────────────
    
    # Join similarity score onto evidence; drop abstract text
    evidence_scored <- merge(
        evidence,
        abstracts_df[, c("pmid", "similarity")],
        by   = "pmid",
        all.x = TRUE
    )
    evidence_scored$similarity[is.na(evidence_scored$similarity)] <- 0
    
    # Keep only evidence rows whose PMID passed
    evidence_filtered <- evidence_scored[
        evidence_scored$pmid %in% passing_pmids,
        c("source", "target", "interaction", "site",
          "evidenceLink", "stmt_hash", "text", "pmid", "similarity")
    ]
    
    # Keep edges that have at least one surviving evidence row
    surviving_hashes <- unique(evidence_filtered$stmt_hash)
    edges_filtered   <- edges[edges$stmt_hash %in% surviving_hashes, ]
    
    # Keep nodes present in surviving edges
    surviving_nodes  <- union(edges_filtered$source, edges_filtered$target)
    nodes_filtered   <- nodes[nodes[[1]] %in% surviving_nodes, ]   # assumes first col is node ID
    
    cat(sprintf(
        "Retained: %d edges (of %d), %d nodes (of %d), %d evidence rows (of %d)\n",
        nrow(edges_filtered),   nrow(edges),
        nrow(nodes_filtered),   nrow(nodes),
        nrow(evidence_filtered), nrow(evidence_scored)
    ))
    
    return(list(
        nodes    = nodes_filtered,
        edges    = edges_filtered,
        evidence = evidence_filtered
    ))
}


# ── Internal helpers ──────────────────────────────────────────────────────────

#' Extract evidence text from edges dataframe via INDRA API
#' @param df Edges dataframe with columns: source, target, interaction, site,
#'           evidenceLink, stmt_hash
#' @return Dataframe with additional columns: text, pmid
#' @keywords internal
#' @noRd
.extract_evidence_text <- function(df) {
    
    required_cols <- c("source", "target", "interaction", "site", "evidenceLink", "stmt_hash")
    missing_cols  <- setdiff(required_cols, names(df))
    if (length(missing_cols) > 0) {
        stop(sprintf("Missing required columns: %s", paste(missing_cols, collapse = ", ")))
    }
    
    results_list <- list()
    result_count <- 0
    unique_hashes <- unique(df$stmt_hash)
    n_hashes      <- length(unique_hashes)
    
    cat(sprintf("Processing %d unique statement hashes...\n", n_hashes))
    
    for (i in seq_along(unique_hashes)) {
        stmt_hash <- unique_hashes[i]
        
        if (i %% 10 == 0) cat(sprintf("Progress: %d/%d\n", i, n_hashes))
        
        evidence_list    <- .query_indra_evidence(stmt_hash)
        if (is.null(evidence_list) || length(evidence_list) == 0) next
        
        matching_indices <- which(df$stmt_hash == stmt_hash)
        
        for (evidence in evidence_list) {
            if (!is.null(evidence[["text"]]) && nchar(evidence[["text"]]) > 0) {
                for (idx in matching_indices) {
                    result_count <- result_count + 1
                    results_list[[result_count]] <- data.frame(
                        source       = df$source[idx],
                        target       = df$target[idx],
                        interaction  = df$interaction[idx],
                        site         = df$site[idx],
                        evidenceLink = df$evidenceLink[idx],
                        stmt_hash    = df$stmt_hash[idx],
                        text         = evidence[["text"]],
                        pmid         = if (is.null(evidence[["pmid"]])) "" else evidence[["pmid"]],
                        stringsAsFactors = FALSE
                    )
                }
            }
        }
    }
    
    if (result_count == 0) {
        warning("No evidence text found for any statement hash")
        return(data.frame(
            source = character(), target = character(), interaction = character(),
            site = character(), evidenceLink = character(), stmt_hash = character(),
            text = character(), pmid = character(), stringsAsFactors = FALSE
        ))
    }
    
    results_df <- do.call(rbind, results_list)
    cat(sprintf("\nComplete! Found %d evidence text entries.\n", nrow(results_df)))
    return(results_df)
}


#' Fetch and clean PubMed abstracts via rentrez
#' @param pmids Character vector of PubMed IDs
#' @return Named list: pmid -> abstract text
#' @keywords internal
#' @importFrom rentrez entrez_fetch
#' @importFrom xml2 read_xml xml_find_first xml_text
#' @noRd
.fetch_clean_abstracts_xml <- function(pmids) {
    results <- list()
    total   <- length(pmids)
    
    cat(sprintf("Fetching %d abstracts...\n", total))
    
    for (i in seq_along(pmids)) {
        pmid <- pmids[i]
        tryCatch({
            record        <- entrez_fetch(db = "pubmed", id = pmid, rettype = "xml")
            doc           <- read_xml(record)
            abstract_node <- xml_find_first(doc, ".//AbstractText")
            
            if (!is.na(abstract_node)) {
                results[[pmid]] <- xml_text(abstract_node)
            }
            
            if (i %% 10 == 0 || i == total) {
                cat(sprintf("Progress: %d/%d (%.1f%%)\n", i, total, (i / total) * 100))
            }
            
            Sys.sleep(0.34)   # respect NCBI rate limit
        }, error = function(e) {
            results[[pmid]] <- ""
            cat(sprintf("Error fetching PMID %s at %d/%d: %s\n", pmid, i, total, e$message))
        })
    }
    
    cat("Done fetching abstracts!\n")
    return(results)
}

#' Query INDRA API for evidence text
#'
#' @param stmt_hash A statement hash string
#' @importFrom httr POST status_code content content_type_json
#' @importFrom jsonlite fromJSON
#' @noRd
#' @return A list of evidence objects from the API, or NULL if error
.query_indra_evidence <- function(stmt_hash) {
    url <- "https://discovery.indra.bio/api/get_evidences_for_stmt_hash"
    
    tryCatch({
        response <- POST(
            url,
            body = list(stmt_hash = stmt_hash),
            encode = "json",
            content_type_json()
        )
        
        if (status_code(response) != 200) {
            warning(sprintf("API returned status %d for stmt_hash: %s", 
                            status_code(response), stmt_hash))
            return(NULL)
        }
        
        content(response, as = "parsed")
    }, error = function(e) {
        warning(sprintf("Error querying stmt_hash %s: %s", stmt_hash, e$message))
        return(NULL)
    })
}