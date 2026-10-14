# Script to create cptac-ccrcc-network.rds in extdata folder
# The PTM-Analysis vignette loads these results instead of calling INDRA
# CoGEx, so that building the vignette needs no network access. Run from
# the package root after changing the vignette's queries.
library(MSstatsBioNet)
ptm = data.table::fread(system.file(
    "extdata/cptac-ccrcc-ptm-adjusted.csv",
    package = "MSstatsBioNet"
))
global = data.table::fread(system.file(
    "extdata/cptac-ccrcc-protein.csv",
    package = "MSstatsBioNet"
))

indra = indra_backend()
sites = convert_ids(indra, prepare_entities(ptm, entity_type = "ptm_site",
                                            id_type = "uniprot"))
proteins = convert_ids(indra, prepare_entities(global,
                                               entity_type = "protein",
                                               id_type = "uniprot"))

# Kinases that act on many substrates, added to the query
kinases = c(MTOR = "HGNC:3942", AKT1 = "HGNC:391", MAPK1 = "HGNC:6871",
            MAPK3 = "HGNC:6877", CDK1 = "HGNC:1722", CSNK2A1 = "HGNC:2457",
            PRKACA = "HGNC:9380", GSK3B = "HGNC:4617", SRC = "HGNC:11283",
            PRKCA = "HGNC:9393")
selected = select_entities(sites, pvalue_cutoff = 0.05)

saved = list(
    sites = sites,
    proteins = proteins,
    network_sites = get_network(indra, selected, subnetwork_query(),
                                interaction_types = "Phosphorylation"),
    network = get_network(indra, selected, subnetwork_query(),
                          interaction_types = "Phosphorylation",
                          include_entities = unname(kinases))
)
saveRDS(saved, "inst/extdata/cptac-ccrcc-network.rds", compress = "xz")
