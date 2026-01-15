test_that("annotateProteinInfoFromIndra works correctly with Uniprot_Mnemonic", {
    df <- data.frame(Protein = c("CLH1_HUMAN"))
    annotated_df <- annotateProteinInfoFromIndra(df, "Uniprot_Mnemonic")
    
    expect_true("Protein" %in% colnames(annotated_df))
    expect_true("UniprotId" %in% colnames(annotated_df))
    expect_true("HgncId" %in% colnames(annotated_df))
    expect_true("HgncName" %in% colnames(annotated_df))
    expect_true("IsTranscriptionFactor" %in% colnames(annotated_df))
    expect_true("IsKinase" %in% colnames(annotated_df))
    expect_true("IsPhosphatase" %in% colnames(annotated_df))

    expect_false(is.na(annotated_df$UniprotId))
    expect_false(is.na(annotated_df$HgncId))
    expect_false(is.na(annotated_df$HgncName))
    expect_false(is.na(annotated_df$IsTranscriptionFactor))
    expect_false(is.na(annotated_df$IsKinase))
    expect_false(is.na(annotated_df$IsPhosphatase))
    
    expect_equal(annotated_df$Protein, "CLH1_HUMAN")
    expect_equal(annotated_df$UniprotId, "Q00610")
    expect_equal(annotated_df$HgncId, "2092")
    expect_equal(annotated_df$HgncName, "CLTC")
    expect_equal(annotated_df$IsTranscriptionFactor, FALSE)
    expect_equal(annotated_df$IsKinase, FALSE)
    expect_equal(annotated_df$IsPhosphatase, FALSE)

})

test_that("annotateProteinInfoFromIndra throws error for missing Protein column", {
    df <- data.frame(NotProtein = c("CLH1_HUMAN"))
    expect_error(annotateProteinInfoFromIndra(df, "Uniprot_Mnemonic"), "Input dataframe must contain 'Protein' column.")
})

test_that("annotateProteinInfoFromIndra works correctly with Uniprot", {
    df <- data.frame(Protein = c("Q00610"))
    annotated_df <- annotateProteinInfoFromIndra(df, "Uniprot")
    
    expect_true("Protein" %in% colnames(annotated_df))
    expect_true("UniprotId" %in% colnames(annotated_df))
    expect_true("HgncId" %in% colnames(annotated_df))
    expect_true("HgncName" %in% colnames(annotated_df))
    expect_true("IsTranscriptionFactor" %in% colnames(annotated_df))
    expect_true("IsKinase" %in% colnames(annotated_df))
    expect_true("IsPhosphatase" %in% colnames(annotated_df))

    expect_false(is.na(annotated_df$UniprotId))
    expect_false(is.na(annotated_df$HgncId))
    expect_false(is.na(annotated_df$HgncName))
    expect_false(is.na(annotated_df$IsTranscriptionFactor))
    expect_false(is.na(annotated_df$IsKinase))
    expect_false(is.na(annotated_df$IsPhosphatase))
    
    expect_equal(annotated_df$Protein, "Q00610")
    expect_equal(annotated_df$UniprotId, "Q00610")
    expect_equal(annotated_df$HgncId, "2092")
    expect_equal(annotated_df$HgncName, "CLTC")
    expect_equal(annotated_df$IsTranscriptionFactor, FALSE)
    expect_equal(annotated_df$IsKinase, FALSE)
    expect_equal(annotated_df$IsPhosphatase, FALSE)
})

test_that("annotateProteinInfoFromIndra returns NA for unknown protein id", {
    df <- data.frame(Protein = c("ABC"))
    annotated_df <- annotateProteinInfoFromIndra(df, "Uniprot_Mnemonic")
    
    expect_true("Protein" %in% colnames(annotated_df))
    expect_true("UniprotId" %in% colnames(annotated_df))
    expect_true("HgncId" %in% colnames(annotated_df))
    expect_true("HgncName" %in% colnames(annotated_df))
    expect_true("IsTranscriptionFactor" %in% colnames(annotated_df))
    expect_true("IsKinase" %in% colnames(annotated_df))
    expect_true("IsPhosphatase" %in% colnames(annotated_df))

    expect_true(is.na(annotated_df$UniprotId))
    expect_true(is.na(annotated_df$HgncId))
    expect_true(is.na(annotated_df$HgncName))
    expect_true(is.na(annotated_df$IsTranscriptionFactor))
    expect_true(is.na(annotated_df$IsKinase))
    expect_true(is.na(annotated_df$IsPhosphatase))
    
    expect_equal(annotated_df$Protein, "ABC")
})

test_that("annotateProteinInfoFromIndra works correctly with HGNC name", {
    df <- data.frame(Protein = c("EGFR"))
    annotated_df <- annotateProteinInfoFromIndra(df, "Hgnc_Name")
    
    expect_true("Protein" %in% colnames(annotated_df))
    expect_true("UniprotId" %in% colnames(annotated_df))
    expect_true("HgncId" %in% colnames(annotated_df))
    expect_true("HgncName" %in% colnames(annotated_df))
    expect_true("IsTranscriptionFactor" %in% colnames(annotated_df))
    expect_true("IsKinase" %in% colnames(annotated_df))
    expect_true("IsPhosphatase" %in% colnames(annotated_df))
    
    expect_true(is.na(annotated_df$UniprotId))
    expect_false(is.na(annotated_df$HgncId))
    expect_false(is.na(annotated_df$HgncName))
    expect_false(is.na(annotated_df$IsTranscriptionFactor))
    expect_false(is.na(annotated_df$IsKinase))
    expect_false(is.na(annotated_df$IsPhosphatase))
    
    expect_equal(annotated_df$Protein, "EGFR")
    expect_equal(annotated_df$HgncId, "3236")
    expect_equal(annotated_df$HgncName, "EGFR")
    expect_type(annotated_df$IsTranscriptionFactor, "logical")
    expect_type(annotated_df$IsKinase, "logical")
    expect_type(annotated_df$IsPhosphatase, "logical")
    
})