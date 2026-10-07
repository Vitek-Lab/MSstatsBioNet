# Script to create panitchpakdi-2022-network.rds in extdata folder
# The Metabolomics-Analysis vignette loads these results instead of calling
# Gilda and INDRA CoGEx, so that building the vignette needs no network
# access. Run from the package root after changing the vignette's queries.
library(MSstatsBioNet)
input = data.table::fread(system.file(
    "extdata/panitchpakdi-2022.csv",
    package = "MSstatsBioNet"
))
input$entity_type = ifelse(input$Analyte == "Diphenhydramine",
                           "drug", "metabolite")
entities = prepare_entities(input,
                            id_column = "Analyte",
                            entity_type = "entity_type",
                            id_type = "chemical_name")

indra = indra_backend()
entities = convert_ids(indra, entities)

proteins = c("HGNC:5182", "HGNC:2625", "HGNC:7967")
selected = select_entities(entities, pvalue_cutoff = 0.05)
selected_inf = select_entities(entities, pvalue_cutoff = 0.05,
                               include_infinite_fc = TRUE)

saved = list(
    entities = entities,
    network_compounds = get_network(indra, selected, subnetwork_query()),
    network = get_network(indra, selected, subnetwork_query(),
                          include_entities = proteins),
    network_inf = get_network(indra, selected_inf, subnetwork_query(),
                              include_entities = proteins)
)
saveRDS(saved, "inst/extdata/panitchpakdi-2022-network.rds")
