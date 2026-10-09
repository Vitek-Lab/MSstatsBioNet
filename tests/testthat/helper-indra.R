# INDRA mocks shared by the backend test files. CoGEx is never called:
# .callIndraCogexApi() returns the saved response in inst/extdata.

.selected_input <- function() {
    input <- data.table::fread(
        system.file("extdata/groupComparisonModel.csv", package = "MSstatsBioNet")
    )
    select_entities(.build_entities_from_annotated_input(input))
}

.mock_indra_response <- function(env = parent.frame()) {
    local_mocked_bindings(
        .callIndraCogexApi = function(ns, ids, fio, cogex_url) {
            readRDS(system.file("extdata/indraResponse.rds", package = "MSstatsBioNet"))
        },
        .get_current_time = function() .fixed_retrieval_time(),
        .env = env
    )
}

.fixed_retrieval_time <- function() {
    as.POSIXct("2026-10-09 12:00:00", tz = "UTC")
}
