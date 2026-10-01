test_that("getSubnetworkFromIndra works correctly", {
    input <- data.table::fread(
        system.file("extdata/groupComparisonModel.csv", package = "MSstatsBioNet")
    )
    local_mocked_bindings(.callIndraCogexApi = function(ns, ids, fio) {
        return(readRDS(system.file("extdata/indraResponse.rds", package = "MSstatsBioNet")))
    })
    suppressWarnings(subnetwork <- getSubnetworkFromIndra(input, statement_types = c("Activation", "Phosphorylation")))
    expect_equal(nrow(subnetwork$nodes), 3)
    expect_equal(nrow(subnetwork$edges), 2)
})

test_that("getSubnetworkFromIndra with different statement type works correctly", {
    input <- data.table::fread(
        system.file("extdata/groupComparisonModel.csv", package = "MSstatsBioNet")
    )
    local_mocked_bindings(.callIndraCogexApi = function(ns, ids, fio) {
        return(readRDS(system.file("extdata/indraResponse.rds", package = "MSstatsBioNet")))
    })
    suppressWarnings(
        subnetwork <- getSubnetworkFromIndra(input, statement_types = c("Complex"))
    )
    expect_equal(nrow(subnetwork$nodes), 8)
    expect_equal(nrow(subnetwork$edges), 16)
})

test_that("Exception is thrown for 400+ proteins in dataframe", {
    input_400 <- data.frame(
        Protein         = paste0("Protein", 1:400),
        log2FC          = rep(1.0, 400),
        adj.pvalue      = rep(0.05, 400),
        EntityNamespace = rep("HGNC", 400),
        EntityId        = paste0("HGNCID", 1:400),
        EntityName      = paste0("HGNCNAME", 1:400),
        issue           = NA
    )
    expect_error(
        getSubnetworkFromIndra(input_400),
        "Invalid Input Error: INDRA query must contain less than 400 proteins.  Consider lowering your p-value cutoff"
    )
})

test_that("Exception is thrown for missing columns in input", {
    input_missing_cols <- data.frame(
        Protein = paste0("Protein", 1:10),
        issue = NA,
        adj.pvalue = 0.05
    )
    expect_error(
        getSubnetworkFromIndra(input_missing_cols),
        "Invalid Input Error: input is missing required column\\(s\\): log2FC, EntityNamespace, EntityId, EntityName\\."
    )
})

# ----- Deprecated arguments (Phase 0 of the API refactor) -----

.run_deprecation_case <- function(...) {
    input <- data.table::fread(
        system.file("extdata/groupComparisonModel.csv", package = "MSstatsBioNet")
    )
    local_mocked_bindings(
        .callIndraCogexApi = function(ns, ids, fio) {
            readRDS(system.file("extdata/indraResponse.rds", package = "MSstatsBioNet"))
        },
        .env = parent.frame()
    )
    getSubnetworkFromIndra(input, statement_types = c("Complex"), ...)
}

test_that("paper_count_cutoff warns as deprecated and is ignored", {
    suppressWarnings(baseline <- .run_deprecation_case())
    suppressWarnings(expect_warning(
        subnetwork <- .run_deprecation_case(paper_count_cutoff = 1),
        "'paper_count_cutoff'.*deprecated"
    ))
    expect_equal(subnetwork, baseline)
    # Values >= 2 used to remove every edge and error; now they are ignored
    suppressWarnings(expect_warning(
        subnetwork <- .run_deprecation_case(paper_count_cutoff = 2),
        "'paper_count_cutoff'.*deprecated"
    ))
    expect_equal(subnetwork, baseline)
})

test_that("correlation_cutoff warns as deprecated", {
    suppressWarnings(expect_warning(
        .run_deprecation_case(correlation_cutoff = 0.5),
        "'correlation_cutoff'.*deprecated"
    ))
})

test_that("protein_level_data warns as deprecated and still adds correlations", {
    input <- data.table::fread(
        system.file("extdata/groupComparisonModel.csv", package = "MSstatsBioNet")
    )
    proteins <- unique(input$Protein)
    protein_level_data <- data.frame(
        Protein        = rep(proteins, each = 3),
        originalRUN    = rep(paste0("run", 1:3), times = length(proteins)),
        LogIntensities = sin(seq_len(3 * length(proteins)))
    )
    suppressWarnings(expect_warning(
        subnetwork <- .run_deprecation_case(
            protein_level_data = protein_level_data, correlation_cutoff = 0
        ),
        "'protein_level_data'.*deprecated"
    ))
    expect_true("correlation" %in% colnames(subnetwork$edges))
})

test_that("default and explicit NULL arguments give no deprecation warning", {
    warns <- capture_warnings(.run_deprecation_case())
    expect_false(any(grepl("deprecated", warns)))
    warns <- capture_warnings(.run_deprecation_case(protein_level_data = NULL))
    expect_false(any(grepl("deprecated", warns)))
})
