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
url = 'https://db.indra.bio/statements/from_agents?subject=7010@HGNC'
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

for (key in keys(edgeToMetadataMapping)) {
    edgeToMetadataMapping[[key]]$data$stmt_type <-
        paste(unique(edgeToMetadataMapping[[key]]$data$stmt_type), 
              collapse = ", ")
}




