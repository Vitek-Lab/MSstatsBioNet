### Fetch Abstracts
library(rentrez)
library(xml2)

fetch_clean_abstracts_xml <- function(pmids) {
    results <- list()
    total <- length(pmids)
    
    # Print initial message
    cat(sprintf("Fetching %d abstracts...\n", total))
    
    for (i in seq_along(pmids)) {
        pmid <- pmids[i]
        
        tryCatch({
            # Fetch as XML (much cleaner to parse)
            record <- entrez_fetch(db = "pubmed", id = pmid, rettype = "xml")
            
            # Parse XML
            doc <- read_xml(record)
            
            # Extract abstract text directly
            abstract_node <- xml_find_first(doc, ".//AbstractText")
            
            if (!is.na(abstract_node)) {
                abstract_text <- xml_text(abstract_node)
                results[[pmid]] <- abstract_text
            }
            
            # Progress update every 10 records or on last record
            if (i %% 10 == 0 || i == total) {
                cat(sprintf("Progress: %d/%d (%.1f%%)\n", i, total, (i/total)*100))
            }
            
            Sys.sleep(0.34)
        }, error = function(e) {
            results[[pmid]] <- paste("Error:", e$message)
            cat(sprintf("Error fetching PMID %s at %d/%d\n", pmid, i, total))
        })
    }
    
    cat("Done!\n")
    return(results)
}

# Usage
input = data.table::fread("~/Downloads/chek1.csv")
annotated_df = annotateProteinInfoFromIndra(input, "Uniprot_Mnemonic")
subnetwork = getSubnetworkFromIndra(
    annotated_df, 
    pvalueCutoff = 0.2,
    logfc_cutoff = NULL,
    evidence_count_cutoff = 1,
    sources_filter = NULL,
    force_include_other = "HGNC:1925",
    filter_by_curation = FALSE
)
pmids <- unique(subnetwork$evidence$pmid)
clean_results <- fetch_clean_abstracts_xml(pmids)
abstracts <- data.table::data.table(
    pmid = names(clean_results),
    abstract = unlist(clean_results, use.names = FALSE)
)

### Get cosine similarity scores between user query and abstracts
library(text2vec)
library(Matrix)

my_query <- "DNA damage repair cancer oncology"

# ── 2. Combine query + abstracts, then vectorise with TF-IDF ─────────────────

all_texts <- c(my_query, abstracts$abstract)

tokens <- itoken(all_texts,
                 preprocessor = tolower,
                 tokenizer    = word_tokenizer)

library(stopwords)
vocab  <- create_vocabulary(tokens, stopwords = stopwords("en"))
vocab  <- prune_vocabulary(vocab, term_count_min = 1)

vectorizer <- vocab_vectorizer(vocab)
dtm        <- create_dtm(tokens, vectorizer)

tfidf      <- TfIdf$new()
dtm_tfidf  <- fit_transform(dtm, tfidf)

# ── 3. Cosine similarity between query (row 1) and each abstract ──────────────

cos_sim <- function(a, b) {
    as.numeric((a %*% t(b)) / (norm(a, "F") * norm(b, "F")))
}

query_vec     <- dtm_tfidf[1, , drop = FALSE]
abstract_vecs <- dtm_tfidf[-1, , drop = FALSE]

scores <- sapply(seq_len(nrow(abstract_vecs)), function(i) {
    cos_sim(query_vec, abstract_vecs[i, , drop = FALSE])
})

# ── 4. Build results data frame ───────────────────────────────────────────────

results <- data.frame(
    abstract        = abstracts,
    similarity      = round(scores, 4),
    stringsAsFactors = FALSE
)

results <- results[order(results$similarity, decreasing = TRUE), ]

# ── 5. Filter by threshold ────────────────────────────────────────────────────

threshold      <- 0.10   # <-- adjust this to taste
filtered       <- results[results$similarity >= threshold, ]

cat("All abstracts with similarity scores:\n")
print(results, row.names = FALSE)

cat(sprintf("\nAbstracts above similarity threshold (%.2f):\n", threshold))
print(filtered, row.names = FALSE)
