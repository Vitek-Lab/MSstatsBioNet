# Test .callGetUniprotIdsFromUniprotMnemonicIdsApi
test_that(".callGetUniprotIdsFromUniprotMnemonicIdsApi works correctly", {
    uniprotMnemonicIds <- list("CLH1_HUMAN")
    local_mocked_bindings(.callGetUniprotIdsFromUniprotMnemonicIdsApi = function(x) {
        return(list(CLH1_HUMAN = "Q00610"))
    })
    result <- .callGetUniprotIdsFromUniprotMnemonicIdsApi(uniprotMnemonicIds)
    expect_type(result, "list")
    expect_true(length(result) == 1)
    expected_value <- list(CLH1_HUMAN = "Q00610")
    expect_equal(result, expected_value)
})

# Test .callGetHgncIdsFromUniprotIdsApi
test_that(".callGetHgncIdsFromUniprotIdsApi works correctly", {
    uniprotIds <- list("Q00610")
    local_mocked_bindings(.callGetHgncIdsFromUniprotIdsApi = function(x) {
        return(list("Q00610" = "2092"))
    })
    result <- .callGetHgncIdsFromUniprotIdsApi(uniprotIds)
    expect_type(result, "list")
    expect_true(length(result) == 1)
    expected_value <- list("Q00610" = "2092")
    expect_equal(result, expected_value)
})

# Test .callGetHgncNamesFromHgncIdsApi
test_that(".callGetHgncNamesFromHgncIdsApi works correctly", {
    hgncIds <- list("2092")
    local_mocked_bindings(.callGetHgncNamesFromHgncIdsApi = function(x) {
        return(list("2092" = "CLTC"))
    })
    result <- .callGetHgncNamesFromHgncIdsApi(hgncIds)
    expect_type(result, "list")
    expect_true(length(result) == 1)
    expected_value <- list("2092" = "CLTC")
    expect_equal(result, expected_value)
})

# Test .callIsKinaseApi
test_that(".callIsKinaseApi works correctly", {
    kinaseGenes <- list("CHEK1")
    local_mocked_bindings(.callIsKinaseApi = function(x) {
        return(list("CHEK1" = TRUE))
    })
    result <- .callIsKinaseApi(kinaseGenes)
    expect_type(result, "list")
    expect_true(length(result) == 1)
    expected_value <- list("CHEK1" = TRUE)
    expect_equal(result, expected_value)
})

# Test .callIsPhosphataseApi
test_that(".callIsPhosphataseApi works correctly", {
    phosphataseGenes <- list("MTM1")
    local_mocked_bindings(.callIsPhosphataseApi = function(x) {
        return(list("MTM1" = TRUE))
    })
    result <- .callIsPhosphataseApi(phosphataseGenes)
    expect_type(result, "list")
    expect_true(length(result) == 1)
    expected_value <- list("MTM1" = TRUE)
    expect_equal(result, expected_value)
})

# Test .callIsTranscriptionFactorApi
test_that(".callIsTranscriptionFactorApi works correctly", {
    transcriptionFactorGenes <- list("STAT1")
    local_mocked_bindings(.callIsTranscriptionFactorApi = function(x) {
        return(list("STAT1" = TRUE))
    })
    result <- .callIsTranscriptionFactorApi(transcriptionFactorGenes)
    expect_type(result, "list")
    expect_true(length(result) == 1)
    expected_value <- list("STAT1" = TRUE)
    expect_equal(result, expected_value)
})

test_that(".callGetHgncIdsFromGildaApi works correctly", {
    hgnc_names <- list("EGFR", "CHEK1")
    result <- .callGetHgncIdsFromGildaApi(hgnc_names)
    expect_type(result, "character")
    expect_true(length(result) == 2)
    expected_value <- c("EGFR" = "3236", "CHEK1" = "1925")
    expect_equal(result, expected_value)
})

