#' Filter a subnetwork by contextual relevance
#'
#' Fetches PubMed abstracts for evidence PMIDs, scores each abstract against a
#' user-supplied query, and returns only the nodes, edges, and evidence rows
#' whose abstracts meet the scoring cutoff.
#'
#' Two scoring methods are available, controlled by the \code{method} argument:
#'
#' \describe{
#'   \item{\code{"tag_count"} (default)}{
#'     Counts how many tags from \code{query} appear as substrings in the
#'     abstract (case-insensitive). The score for each abstract is an integer
#'     in \code{[0, length(query)]}. Set \code{cutoff} to the minimum number of
#'     tags that must appear - e.g. \code{cutoff = 2} keeps abstracts that
#'     mention at least 2 of your tags. \code{query} must be a character
#'     \emph{vector} of tags when using this method.
#'   }
#'   \item{\code{"cosine"}}{
#'     Scores abstracts using TF-IDF cosine similarity against \code{query}.
#'     Scores are in \code{[-1, 1]} (in practice \code{[0, 1]} for text).
#'     Set \code{cutoff} to a decimal threshold - e.g. \code{cutoff = 0.10}.
#'     \code{query} should be a single character string; expand it with
#'     synonyms and related terms for better recall under exact token matching.
#'   }
#' }
#'
#' @param nodes  A dataframe of network nodes.
#' @param edges  A dataframe of network edges with columns: source, target,
#'               interaction, site, evidenceLink, stmt_hash.
#' @param query  For \code{method = "tag_count"}: a character vector of tags,
#'               e.g. \code{c("CHEK1", "DNA damage", "DNA damage repair")}.
#'               For \code{method = "cosine"}: a single character string.
#' @param cutoff Numeric threshold applied to the chosen scoring method.
#'               \itemize{
#'                 \item \code{"tag_count"}: integer >= 0; abstracts must
#'                   contain at least this many tags. Max possible value is
#'                   \code{length(query)}. Default \code{1}.
#'                 \item \code{"cosine"}: numeric in \code{[-1, 1]}; abstracts
#'                   must score >= this value. Default \code{0.10}.
#'               }
#' @param method One of \code{"tag_count"} (default) or \code{"cosine"}.
#'
#' @return A named list with three elements:
#'   \item{nodes}{Filtered nodes dataframe (only nodes present in kept edges)}
#'   \item{edges}{Filtered edges dataframe}
#'   \item{evidence}{Dataframe with columns: source, target, interaction, site,
#'     evidenceLink, stmt_hash, text, pmid, score. The \code{score} column
#'     contains tag counts (integer) or cosine similarities (numeric) depending
#'     on the method used.}
#'
#' @importFrom text2vec itoken word_tokenizer create_vocabulary prune_vocabulary vocab_vectorizer create_dtm TfIdf fit_transform
#' @importFrom stopwords stopwords
#' @note **Beta feature:** This function is experimental and the API may
#'   change without notice in future versions.
#' @export
filterSubnetworkByContext <- function(nodes,
                                      edges,
                                      query,
                                      cutoff = NULL,
                                      method = c("tag_count", "cosine")) {
    
    method <- match.arg(method)

    if (method == "tag_count") {
        if (!is.character(query) || length(query) < 1 || 
            any(is.na(query)) || any(!nzchar(query))) {
            stop("`query` must be a character vector of tags when method = 'tag_count'.")
        }
        if (is.null(cutoff)) cutoff <- 1L
        if (!is.numeric(cutoff) || length(cutoff) != 1L || is.na(cutoff) ||
            cutoff < 0 || cutoff > length(query) || cutoff != as.integer(cutoff)) {
            stop("`cutoff` must be a single integer in [0, length(query)] when method = 'tag_count'.")
        }
        cat(sprintf(
            "Method: tag_count | Tags: %d | Cutoff: >= %d tag(s)\n",
            length(query), cutoff
        ))
    } else {
        if (!is.character(query) || length(query) != 1L ||
            is.na(query) || !nzchar(query)) {
            stop("`query` must be a single character string when method = 'cosine'.")
        }
        if (is.null(cutoff)) cutoff <- 0.10
        if (!is.numeric(cutoff) || length(cutoff) != 1L || is.na(cutoff) ||
            !is.finite(cutoff) || cutoff < 0 || cutoff > 1) {
            stop("`cutoff` must be a single numeric value in [0, 1] when method = 'cosine'.")
        }
        cat(sprintf(
            "Method: cosine | Cutoff: >= %.2f\n", cutoff
        ))
    }
    
    evidence <- .extract_evidence_text(edges)
    
    if (nrow(evidence) == 0) {
        evidence$score <- if (method == "tag_count") integer(0) else numeric(0)
        warning("No evidence text found - returning unfiltered inputs.")
        return(list(nodes = nodes, edges = edges, evidence = evidence))
    }
    pmids <- unique(evidence$pmid[!is.na(evidence$pmid) & nchar(evidence$pmid) > 0])
    
    if (length(pmids) == 0) {
        evidence$score <- if (method == "tag_count") {
            rep(NA_integer_, nrow(evidence))
        } else {
            rep(NA_real_, nrow(evidence))
        }
        warning("No PMIDs found in evidence - returning unfiltered inputs.")
        return(list(nodes = nodes, edges = edges, evidence = evidence))
    }
    
    abstract_list <- .fetch_clean_abstracts_xml(pmids)
    abstracts_df  <- data.frame(
        pmid     = names(abstract_list),
        abstract = unlist(abstract_list, use.names = FALSE),
        stringsAsFactors = FALSE
    )
    
    if (method == "tag_count") {
        abstracts_df$score <- .score_by_tag_count(abstracts_df$abstract, query)
    } else {
        abstracts_df$score <- .score_by_cosine(query, abstracts_df$abstract)
    }
    
    passing_pmids <- abstracts_df$pmid[abstracts_df$score >= cutoff]
    
    cat(sprintf(
        "\n%d / %d abstracts passed cutoff (score >= %s)\n",
        length(passing_pmids), nrow(abstracts_df), cutoff
    ))
    
    evidence_scored <- merge(
        evidence,
        abstracts_df[, c("pmid", "score")],
        by    = "pmid",
        all.x = TRUE
    )
    evidence_scored$score[is.na(evidence_scored$score)] <- 0
    
    evidence_filtered <- evidence_scored[
        evidence_scored$pmid %in% passing_pmids,
        c("source", "target", "interaction", "site",
          "evidenceLink", "stmt_hash", "text", "pmid", "score")
    ]
    
    surviving_hashes <- unique(evidence_filtered$stmt_hash)
    edges_filtered   <- edges[edges$stmt_hash %in% surviving_hashes, ]
    
    surviving_nodes  <- union(edges_filtered$source, edges_filtered$target)
    if (!"id" %in% names(nodes)) {
        stop("`nodes` must contain an `id` column.")
    }
    nodes_filtered   <- nodes[nodes$id %in% surviving_nodes, ]
    
    cat(sprintf(
        "Retained: %d edges (of %d), %d nodes (of %d), %d evidence rows (of %d)\n",
        nrow(edges_filtered),    nrow(edges),
        nrow(nodes_filtered),    nrow(nodes),
        nrow(evidence_filtered), nrow(evidence_scored)
    ))
    
    return(list(
        nodes    = nodes_filtered,
        edges    = edges_filtered,
        evidence = evidence_filtered
    ))
}


#' Score abstracts by tag count
#'
#' For each abstract, counts how many tags appear as case-insensitive substrings.
#'
#' @param abstracts Character vector of abstract texts.
#' @param tags      Character vector of tags to search for.
#' @return Integer vector of tag hit counts, same length as \code{abstracts}.
#' @keywords internal
#' @noRd
.score_by_tag_count <- function(abstracts, tags) {
    abstracts_lower <- tolower(abstracts)
    tags_lower      <- tolower(tags)
    
    sapply(abstracts_lower, function(abstract) {
        sum(sapply(tags_lower, function(tag) grepl(tag, abstract, fixed = TRUE)))
    }, USE.NAMES = FALSE)
}


#' Score abstracts by TF-IDF cosine similarity
#'
#' Vectorises the query and all abstracts together with TF-IDF, then returns
#' the cosine similarity of each abstract against the query vector.
#'
#' @param query     Single character string query.
#' @param abstracts Character vector of abstract texts.
#' @return Numeric vector of cosine similarities in \code{[0, 1]}, same length
#'   as \code{abstracts}.
#' @keywords internal
#' @noRd
#' @importFrom text2vec itoken word_tokenizer create_vocabulary prune_vocabulary vocab_vectorizer create_dtm TfIdf fit_transform
#' @importFrom stopwords stopwords
.score_by_cosine <- function(query, abstracts) {
    all_texts  <- c(query, abstracts)
    
    tokens     <- itoken(all_texts,
                         preprocessor = tolower,
                         tokenizer    = word_tokenizer)
    vocab      <- create_vocabulary(tokens, stopwords = stopwords("en"))
    vocab      <- prune_vocabulary(vocab, term_count_min = 1)
    vectorizer <- vocab_vectorizer(vocab)
    dtm        <- create_dtm(tokens, vectorizer)
    tfidf      <- TfIdf$new()
    dtm_tfidf  <- fit_transform(dtm, tfidf)
    
    .cos_sim <- function(a, b) {
        a     <- as.numeric(a)
        b     <- as.numeric(b)
        denom <- sqrt(sum(a^2)) * sqrt(sum(b^2))
        if (denom == 0) return(0)
        sum(a * b) / denom
    }
    
    query_vec     <- dtm_tfidf[1, , drop = FALSE]
    abstract_vecs <- dtm_tfidf[-1, , drop = FALSE]
    
    scores <- sapply(seq_len(nrow(abstract_vecs)), function(i) {
        .cos_sim(query_vec, abstract_vecs[i, , drop = FALSE])
    })
    
    round(scores, 4)
}



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
    
    results_list  <- list()
    result_count  <- 0
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
#' @noRd
#' @importFrom rentrez entrez_fetch
#' @importFrom xml2 read_xml xml_find_all xml_text
.fetch_clean_abstracts_xml <- function(pmids) {
    results <- list()
    total   <- length(pmids)
    
    cat(sprintf("Fetching %d abstracts...\n", total))
    
    for (i in seq_along(pmids)) {
        pmid <- pmids[i]
        
        record <- tryCatch(
            entrez_fetch(db = "pubmed", id = pmid, rettype = "xml"),
            error = function(e) {
                cat(sprintf("Error fetching PMID %s at %d/%d: %s\n", pmid, i, total, e$message))
                NULL
            }
        )
        
        if (is.null(record)) {
            results[[pmid]] <- ""
            next
        }
        
        doc <- read_xml(record)
        abstract_nodes <- xml_find_all(doc, ".//AbstractText")
        
        if (length(abstract_nodes) > 0) {
            results[[pmid]] <- paste(trimws(xml_text(abstract_nodes)), collapse = " ")
        } else {
            results[[pmid]] <- ""
        }
        
        if (i %% 10 == 0 || i == total) {
            cat(sprintf("Progress: %d/%d (%.1f%%)\n", i, total, (i / total) * 100))
        }
        
        Sys.sleep(0.34)
    }
    
    cat("Done fetching abstracts!\n")
    return(results)
}


#' Query INDRA API for evidence text
#' @param stmt_hash A statement hash string
#' @return A list of evidence objects from the API, or NULL if error
#' @keywords internal
#' @noRd
#' @importFrom httr POST status_code content content_type_json
#' @importFrom jsonlite fromJSON
.query_indra_evidence <- function(stmt_hash) {
    url <- "https://discovery.indra.bio/api/get_evidences_for_stmt_hash"
    
    tryCatch({
        response <- POST(
            url,
            body   = list(stmt_hash = stmt_hash),
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