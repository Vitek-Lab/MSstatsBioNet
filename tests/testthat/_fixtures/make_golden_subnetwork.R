# Regenerates golden_subnetwork.rds, the getSubnetworkFromIndra() output
# pinned after Phase 1 of the API refactor. Later phases must reproduce it,
# so regenerate only when a change to the output is intended.
#
# Run from the package root:
#   Rscript tests/testthat/_fixtures/make_golden_subnetwork.R
devtools::load_all(quiet = TRUE)
input <- data.table::fread(
    system.file("extdata/groupComparisonModel.csv", package = "MSstatsBioNet")
)
testthat::local_mocked_bindings(
    .callIndraCogexApi = function(ns, ids, fio, cogex_url) {
        readRDS(system.file("extdata/indraResponse.rds", package = "MSstatsBioNet"))
    }
)
subnetwork <- suppressWarnings(getSubnetworkFromIndra(input))
saveRDS(subnetwork, "tests/testthat/_fixtures/golden_subnetwork.rds")
