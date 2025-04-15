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