test_that("getSubnetworkFromIndra works correctly", {
    input <- data.table::fread(
        system.file("extdata/groupComparisonModel.csv", package = "MSstatsBioNet")
    )
    local_mocked_bindings(.callIndraCogexApi = function(x,y) {
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
    local_mocked_bindings(.callIndraCogexApi = function(x,y) {
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
        Protein = paste0("Protein", 1:400),
        HgncId = paste0("HGNCID", 1:400),
        issue = NA,
        adj.pvalue = 0.05
    )
    expect_error(
        getSubnetworkFromIndra(input_400),
        "Invalid Input Error: INDRA query must contain less than 400 proteins.  Consider lowering your p-value cutoff"
    )
})

test_that("Exception is thrown for missing columns in input", {
    input_missing_hgnc_id <- data.frame(
        Protein = paste0("Protein", 1:10),
        issue = NA,
        adj.pvalue = 0.05
    )
    expect_error(
        getSubnetworkFromIndra(input_missing_hgnc_id),
        "Invalid Input Error: Input must contain a column named 'HgncId'."
    )
})
