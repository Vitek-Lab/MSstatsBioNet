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

test_that(".callGroundEntitiesFromGildaApi returns aligned (ns, id, name) per input (live)", {
    text_inputs <- list("EGFR", "CHEK1")
    result <- .callGroundEntitiesFromGildaApi(text_inputs, keep_only = "HGNC")
    expect_type(result, "list")
    expect_true(length(result) == 2)
    expect_setequal(names(result), c("EGFR", "CHEK1"))
    for (input_text in names(result)) {
        g <- result[[input_text]]
        expect_true(all(c("ns", "id", "name") %in% names(g)))
        expect_equal(length(g$ns), length(g$id))
        expect_equal(length(g$ns), length(g$name))
        expect_true(all(g$ns == "HGNC"))
    }
    expect_true("3236"  %in% result[["EGFR"]]$id)
    expect_true("1925"  %in% result[["CHEK1"]]$id)
})

test_that(".callGroundEntitiesFromGildaApi keeps non-HGNC namespaces when keep_only is NULL (mocked)", {
    text_inputs <- list("EGFR", "glucose")
    local_mocked_bindings(.callGroundEntitiesFromGildaApi = function(textInputs, keep_only = NULL) {
        list(
            EGFR    = list(ns = "HGNC",         id = "3236",  name = "EGFR"),
            glucose = list(ns = c("MESH", "CHEBI"),
                           id = c("3815", "17234"),
                           name = c("KIT",  "glucose"))
        )
    })
    result <- .callGroundEntitiesFromGildaApi(text_inputs)
    expect_setequal(names(result), c("EGFR", "glucose"))
    expect_equal(result[["EGFR"]]$ns, "HGNC")
    expect_equal(result[["glucose"]]$ns, c("MESH", "CHEBI"))
    expect_equal(result[["glucose"]]$id, c("3815", "17234"))
    expect_equal(result[["glucose"]]$name, c("KIT", "glucose"))
})

