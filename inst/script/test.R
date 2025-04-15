library(tidyverse)
library(MSstatsBioNet)
input <- data.table::fread(system.file(
    "extdata/model.csv",
    package = "MSstatsBioNet"
))

input = input %>% filter(Label == "DMSO-VTP50469") %>% filter(is.na(issue))
annotated_df = annotateProteinInfoFromIndra(input, "Uniprot_Mnemonic")
main_target = "MEN1_HUMAN"



# Fit a normal distribution to log2FC column values
library(MASS)
log2fc_values <- annotated_df$log2FC
fit <- fitdistr(log2fc_values, "normal")
para <- fit$estimate
pnorm(0, mean = para[1], sd = para[2])

# Fit a negative binomial distribution with parameters
n = 0.5
p = 0.06
probability <- dnbinom(2, size = n, prob = p)
probability

# Call INDRA
library(httr)
url = 'https://db.indra.bio/statements/from_agents?subject=7010@HGNC&max_stmts=1000'
response <- GET(url)
z = content(response)

library(r2r)
edgeToMetadataMapping <- hashmap()

for (index in seq(1, length(z$statements))) {
    edge <- z$statements[[index]]
    key <- ""
    if (edge$type == "Complex") {
        next
    } else if (!("HGNC" %in% names(edge$obj$db_refs))) {
        next
    } else if (!(edge$obj$db_refs$HGNC %in% annotated_df$HgncId)) {
        next
    } else {
        key <- paste(edge$subj$db_refs$HGNC, edge$obj$db_refs$HGNC, sep = "_")
    }
    
    if (key %in% keys(edgeToMetadataMapping)) {
        edgeToMetadataMapping[[key]]$data$evidence_count <-
            edgeToMetadataMapping[[key]]$data$evidence_count +
            z$evidence_counts[[index]]
        edgeToMetadataMapping[[key]]$data$stmt_type <- unique(c(
            edgeToMetadataMapping[[key]]$data$stmt_type,
            edge$type))
    } else {
        # edge <- MSstatsBioNet:::.addAdditionalMetadataToIndraEdge(edge, annotated_df)
        edgeToMetadataMapping[[key]] <- edge
        edgeToMetadataMapping[[key]]$data$evidence_count <-
            z$evidence_counts[[index]]
        edgeToMetadataMapping[[key]]$data$stmt_type <- c(edge$type)
        edgeToMetadataMapping[[key]]$source_id <- edge$subj$db_refs$HGNC
        edgeToMetadataMapping[[key]]$target_id <- edge$obj$db_refs$HGNC
    }
}

# Calculate probabilities
for (key in keys(edgeToMetadataMapping)) {
    edgeToMetadataMapping[[key]]$data$stmt_type <-
        paste(unique(edgeToMetadataMapping[[key]]$data$stmt_type), 
              collapse = ", ")
    prob_logFC = 0
    logFC = annotated_df %>% filter(HgncId == edgeToMetadataMapping[[key]]$target_id)
    logFC = logFC$log2FC[[1]]
    if (logFC > para[1]) {
        prob_logFC = 1 - pnorm(logFC, mean = para[1], sd = para[2])
    } else {
        prob_logFC = pnorm(logFC, mean = para[1], sd = para[2])
    }
    evidence_prob = dnbinom(min(10, edgeToMetadataMapping[[key]]$data$evidence_count), size = n, prob = p)
    edgeToMetadataMapping[[key]]$data$total_prob = prob_logFC * evidence_prob
    edgeToMetadataMapping[[key]]$data$logFC = logFC
}

# Construct DF and sort
edges <- data.frame(
    source = vapply(keys(edgeToMetadataMapping), function(x) {
        query(edgeToMetadataMapping, x)$source_id
    }, ""),
    target = vapply(keys(edgeToMetadataMapping), function(x) {
        query(edgeToMetadataMapping, x)$target_id
    }, ""),
    interaction = vapply(keys(edgeToMetadataMapping), function(x) {
        query(edgeToMetadataMapping, x)$data$stmt_type
    }, ""),
    evidenceCount = vapply(keys(edgeToMetadataMapping), function(x) {
        query(edgeToMetadataMapping, x)$data$evidence_count
    }, 1),
    logFC = vapply(keys(edgeToMetadataMapping), function(x) {
        query(edgeToMetadataMapping, x)$data$logFC
    }, 1),
    prob = vapply(keys(edgeToMetadataMapping), function(x) {
        query(edgeToMetadataMapping, x)$data$total_prob
    }, 1),
    stringsAsFactors = FALSE
)


