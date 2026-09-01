test_that("annotateProteinInfoFromIndra works correctly with Uniprot_Mnemonic", {
    df <- data.frame(Protein = c("CLH1_HUMAN"))
    annotated_df <- annotateProteinInfoFromIndra(df, "Uniprot_Mnemonic")

    expect_true("Protein" %in% colnames(annotated_df))
    expect_true("UniprotId" %in% colnames(annotated_df))
    expect_true("EntityNamespace" %in% colnames(annotated_df))
    expect_true("EntityId" %in% colnames(annotated_df))
    expect_true("EntityName" %in% colnames(annotated_df))
    expect_true("IsTranscriptionFactor" %in% colnames(annotated_df))
    expect_true("IsKinase" %in% colnames(annotated_df))
    expect_true("IsPhosphatase" %in% colnames(annotated_df))

    expect_false(is.na(annotated_df$UniprotId))
    expect_false(is.na(annotated_df$EntityNamespace))
    expect_false(is.na(annotated_df$EntityId))
    expect_false(is.na(annotated_df$EntityName))
    expect_false(is.na(annotated_df$IsTranscriptionFactor))
    expect_false(is.na(annotated_df$IsKinase))
    expect_false(is.na(annotated_df$IsPhosphatase))

    expect_equal(annotated_df$Protein, "CLH1_HUMAN")
    expect_equal(annotated_df$UniprotId, "Q00610")
    expect_equal(annotated_df$EntityNamespace, "HGNC")
    expect_equal(annotated_df$EntityId, "2092")
    expect_equal(annotated_df$EntityName, "CLTC")
    expect_equal(annotated_df$IsTranscriptionFactor, FALSE)
    expect_equal(annotated_df$IsKinase, FALSE)
    expect_equal(annotated_df$IsPhosphatase, FALSE)

})

test_that("annotateProteinInfoFromIndra throws error for missing Protein column", {
    df <- data.frame(NotProtein = c("CLH1_HUMAN"))
    expect_error(annotateProteinInfoFromIndra(df, "Uniprot_Mnemonic"), "Input dataframe must contain 'Protein' column.")
})

test_that("annotateProteinInfoFromIndra throws error for invalid proteinIdType", {
    df <- data.frame(Protein = c("CLH1_HUMAN"))
    expect_error(
        annotateProteinInfoFromIndra(df, "NotAType"),
        "Invalid proteinIdType"
    )
})

test_that("annotateProteinInfoFromIndra works correctly with Uniprot", {
    df <- data.frame(Protein = c("Q00610"))
    annotated_df <- annotateProteinInfoFromIndra(df, "Uniprot")

    expect_true("Protein" %in% colnames(annotated_df))
    expect_true("UniprotId" %in% colnames(annotated_df))
    expect_true("EntityNamespace" %in% colnames(annotated_df))
    expect_true("EntityId" %in% colnames(annotated_df))
    expect_true("EntityName" %in% colnames(annotated_df))
    expect_true("IsTranscriptionFactor" %in% colnames(annotated_df))
    expect_true("IsKinase" %in% colnames(annotated_df))
    expect_true("IsPhosphatase" %in% colnames(annotated_df))

    expect_false(is.na(annotated_df$UniprotId))
    expect_false(is.na(annotated_df$EntityNamespace))
    expect_false(is.na(annotated_df$EntityId))
    expect_false(is.na(annotated_df$EntityName))
    expect_false(is.na(annotated_df$IsTranscriptionFactor))
    expect_false(is.na(annotated_df$IsKinase))
    expect_false(is.na(annotated_df$IsPhosphatase))

    expect_equal(annotated_df$Protein, "Q00610")
    expect_equal(annotated_df$UniprotId, "Q00610")
    expect_equal(annotated_df$EntityNamespace, "HGNC")
    expect_equal(annotated_df$EntityId, "2092")
    expect_equal(annotated_df$EntityName, "CLTC")
    expect_equal(annotated_df$IsTranscriptionFactor, FALSE)
    expect_equal(annotated_df$IsKinase, FALSE)
    expect_equal(annotated_df$IsPhosphatase, FALSE)
})

test_that("annotateProteinInfoFromIndra returns NA for unknown protein id", {
    df <- data.frame(Protein = c("ABC"))
    annotated_df <- annotateProteinInfoFromIndra(df, "Uniprot_Mnemonic")

    expect_true("Protein" %in% colnames(annotated_df))
    expect_true("UniprotId" %in% colnames(annotated_df))
    expect_true("EntityNamespace" %in% colnames(annotated_df))
    expect_true("EntityId" %in% colnames(annotated_df))
    expect_true("EntityName" %in% colnames(annotated_df))
    expect_true("IsTranscriptionFactor" %in% colnames(annotated_df))
    expect_true("IsKinase" %in% colnames(annotated_df))
    expect_true("IsPhosphatase" %in% colnames(annotated_df))

    expect_true(is.na(annotated_df$UniprotId))
    expect_true(is.na(annotated_df$EntityNamespace))
    expect_true(is.na(annotated_df$EntityId))
    expect_true(is.na(annotated_df$EntityName))
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
    expect_true("EntityNamespace" %in% colnames(annotated_df))
    expect_true("EntityId" %in% colnames(annotated_df))
    expect_true("EntityName" %in% colnames(annotated_df))
    expect_true("IsTranscriptionFactor" %in% colnames(annotated_df))
    expect_true("IsKinase" %in% colnames(annotated_df))
    expect_true("IsPhosphatase" %in% colnames(annotated_df))

    expect_true(is.na(annotated_df$UniprotId))
    expect_false(is.na(annotated_df$EntityNamespace))
    expect_false(is.na(annotated_df$EntityId))
    expect_false(is.na(annotated_df$EntityName))
    expect_false(is.na(annotated_df$IsTranscriptionFactor))
    expect_false(is.na(annotated_df$IsKinase))
    expect_false(is.na(annotated_df$IsPhosphatase))

    expect_equal(annotated_df$Protein, "EGFR")
    expect_equal(annotated_df$EntityNamespace, "HGNC")
    expect_equal(annotated_df$EntityId, "3236")
    expect_equal(annotated_df$EntityName, "EGFR")
    expect_type(annotated_df$IsTranscriptionFactor, "logical")
    expect_type(annotated_df$IsKinase, "logical")
    expect_type(annotated_df$IsPhosphatase, "logical")

})

# ----- Protein groups: ";"-joined identifiers are grounded member by member,
# then pooled into the semicolon-joined Entity* columns that
# getSubnetworkFromIndra already fans out. Mocked, so no network access. -----

test_that(".splitProteinGroup splits, trims and drops empty members", {
    expect_equal(MSstatsBioNet:::.splitProteinGroup("P13747;P23132"),
                 c("P13747", "P23132"))
    expect_equal(MSstatsBioNet:::.splitProteinGroup("P13747"), "P13747")
    expect_equal(MSstatsBioNet:::.splitProteinGroup(" P13747 ; P23132 "),
                 c("P13747", "P23132"))
    expect_equal(MSstatsBioNet:::.splitProteinGroup(";P13747;;"), "P13747")
    expect_equal(MSstatsBioNet:::.splitProteinGroup(NA), character(0))
    expect_equal(MSstatsBioNet:::.splitProteinGroup(""), character(0))
})

test_that(".joinProteinGroup round-trips and returns NA when empty", {
    expect_equal(MSstatsBioNet:::.joinProteinGroup(c("P13747", "P23132")),
                 "P13747;P23132")
    expect_equal(MSstatsBioNet:::.joinProteinGroup("P13747"), "P13747")
    expect_true(is.na(MSstatsBioNet:::.joinProteinGroup(character(0))))
})

test_that(".stripPtmSite removes only the site suffix", {
    expect_equal(MSstatsBioNet:::.stripPtmSite(c("P13747_S148", "P23132")),
                 c("P13747", "P23132"))
    expect_equal(MSstatsBioNet:::.stripPtmSite("CLH1_HUMAN"), "CLH1_HUMAN")
})

test_that("PTM site suffixes are stripped from every protein group member", {
    df <- data.frame(Protein = c("P13747_S148;P23132_T20", "P13747_S148"),
                     stringsAsFactors = FALSE)
    out <- MSstatsBioNet:::.populateUniprotIdsInDataFrame(df, "Uniprot")

    expect_equal(out$GlobalProtein, c("P13747;P23132", "P13747"))
    expect_equal(out$UniprotId,     c("P13747;P23132", "P13747"))
})

test_that("annotateProteinInfoFromIndra grounds each Uniprot protein group member", {
    df <- data.frame(Protein = c("P13747;P23132", "Q00610"),
                     stringsAsFactors = FALSE)
    local_mocked_bindings(
        .callGetHgncIdsFromUniprotIdsApi = function(uniprotIds) {
            list(P13747 = "4931", P23132 = "10012", Q00610 = "2092")
        },
        .callGetHgncNamesFromHgncIdsApi = function(hgncIds) {
            list(`4931` = "HLA-E", `10012` = "RAD23A", `2092` = "CLTC")
        },
        .callIsTranscriptionFactorApi = function(genes) list(CLTC = FALSE),
        .callIsKinaseApi = function(genes) list(CLTC = FALSE),
        .callIsPhosphataseApi = function(genes) list(CLTC = FALSE)
    )
    annotated_df <- annotateProteinInfoFromIndra(df, "Uniprot")

    group_row <- annotated_df[annotated_df$Protein == "P13747;P23132", ]
    expect_equal(group_row$UniprotId,       "P13747;P23132")
    expect_equal(group_row$EntityNamespace, "HGNC;HGNC")
    expect_equal(group_row$EntityId,        "4931;10012")
    expect_equal(group_row$EntityName,      "HLA-E;RAD23A")
    # Gene-only flags stay NA while the row carries more than one grounding
    expect_true(is.na(group_row$IsTranscriptionFactor))
    expect_true(is.na(group_row$IsKinase))
    expect_true(is.na(group_row$IsPhosphatase))

    # A single-identifier row is unaffected by the split
    single_row <- annotated_df[annotated_df$Protein == "Q00610", ]
    expect_equal(single_row$EntityNamespace, "HGNC")
    expect_equal(single_row$EntityId,        "2092")
    expect_equal(single_row$EntityName,      "CLTC")
    expect_equal(single_row$IsKinase,        FALSE)
})

test_that("a protein group whose members share a gene collapses to one grounding", {
    df <- data.frame(Protein = "P13747;P13747-2", stringsAsFactors = FALSE)
    local_mocked_bindings(
        .callGetHgncIdsFromUniprotIdsApi = function(uniprotIds) {
            list(P13747 = "4931", `P13747-2` = "4931")
        },
        .callGetHgncNamesFromHgncIdsApi = function(hgncIds) list(`4931` = "HLA-E"),
        .callIsTranscriptionFactorApi = function(genes) list(`HLA-E` = FALSE),
        .callIsKinaseApi = function(genes) list(`HLA-E` = FALSE),
        .callIsPhosphataseApi = function(genes) list(`HLA-E` = TRUE)
    )
    annotated_df <- annotateProteinInfoFromIndra(df, "Uniprot")

    expect_equal(annotated_df$EntityNamespace, "HGNC")
    expect_equal(annotated_df$EntityId,        "4931")
    expect_equal(annotated_df$EntityName,      "HLA-E")
    # Collapsed to a single grounding, so the gene-only flags are populated
    expect_equal(annotated_df$IsTranscriptionFactor, FALSE)
    expect_equal(annotated_df$IsPhosphatase,         TRUE)
})

test_that("unresolvable protein group members are dropped, not carried as NA", {
    df <- data.frame(Protein = c("P13747;NOTANID", "NOTANID;ALSONOT"),
                     stringsAsFactors = FALSE)
    local_mocked_bindings(
        .callGetHgncIdsFromUniprotIdsApi = function(uniprotIds) list(P13747 = "4931"),
        .callGetHgncNamesFromHgncIdsApi = function(hgncIds) list(`4931` = "HLA-E"),
        .callIsTranscriptionFactorApi = function(genes) list(`HLA-E` = FALSE),
        .callIsKinaseApi = function(genes) list(`HLA-E` = FALSE),
        .callIsPhosphataseApi = function(genes) list(`HLA-E` = FALSE)
    )
    annotated_df <- annotateProteinInfoFromIndra(df, "Uniprot")

    # Partially resolved group keeps only the member that grounded, and so is
    # single-grounded and does get the gene-only flags
    expect_equal(annotated_df$EntityId[1],   "4931")
    expect_equal(annotated_df$EntityName[1], "HLA-E")
    expect_equal(annotated_df$IsKinase[1],   FALSE)

    # Fully unresolved group stays NA across the Entity columns, and its
    # presence alongside a resolved row must not break the flag lookups
    expect_true(is.na(annotated_df$EntityNamespace[2]))
    expect_true(is.na(annotated_df$EntityId[2]))
    expect_true(is.na(annotated_df$EntityName[2]))
    expect_true(is.na(annotated_df$IsKinase[2]))
})

test_that("Uniprot_Mnemonic groups map each member to its own UniProt id", {
    df <- data.frame(Protein = c("CLH1_HUMAN;HLAE_HUMAN", "CLH1_HUMAN;NOPE_HUMAN"),
                     stringsAsFactors = FALSE)
    local_mocked_bindings(
        .callGetUniprotIdsFromUniprotMnemonicIdsApi = function(uniprotMnemonicIds) {
            list(CLH1_HUMAN = "Q00610", HLAE_HUMAN = "P13747")
        },
        .callGetHgncIdsFromUniprotIdsApi = function(uniprotIds) {
            list(Q00610 = "2092", P13747 = "4931")
        },
        .callGetHgncNamesFromHgncIdsApi = function(hgncIds) {
            list(`2092` = "CLTC", `4931` = "HLA-E")
        },
        .callIsTranscriptionFactorApi = function(genes) list(CLTC = FALSE),
        .callIsKinaseApi = function(genes) list(CLTC = FALSE),
        .callIsPhosphataseApi = function(genes) list(CLTC = FALSE)
    )
    annotated_df <- annotateProteinInfoFromIndra(df, "Uniprot_Mnemonic")

    expect_equal(annotated_df$UniprotId[1],       "Q00610;P13747")
    expect_equal(annotated_df$EntityNamespace[1], "HGNC;HGNC")
    expect_equal(annotated_df$EntityId[1],        "2092;4931")
    expect_equal(annotated_df$EntityName[1],      "CLTC;HLA-E")

    # Only one member resolves, so the row degrades to that single grounding
    expect_equal(annotated_df$UniprotId[2],  "Q00610")
    expect_equal(annotated_df$EntityName[2], "CLTC")
    expect_equal(annotated_df$IsKinase[2],   FALSE)
})

test_that("Hgnc_Name groups pool and deduplicate Gilda groundings", {
    df <- data.frame(Protein = c("EGFR;ERBB2", "EGFR;EGFR"),
                     stringsAsFactors = FALSE)
    local_mocked_bindings(
        .callGroundEntitiesFromGildaApi = function(textInputs, keep_only = NULL, organisms = NULL) {
            list(EGFR  = list(ns = "HGNC", id = "3236", name = "EGFR"),
                 ERBB2 = list(ns = "HGNC", id = "3430", name = "ERBB2"))
        },
        .callIsTranscriptionFactorApi = function(genes) list(EGFR = FALSE),
        .callIsKinaseApi = function(genes) list(EGFR = TRUE),
        .callIsPhosphataseApi = function(genes) list(EGFR = FALSE)
    )
    annotated_df <- annotateProteinInfoFromIndra(df, "Hgnc_Name")

    expect_equal(annotated_df$EntityNamespace[1], "HGNC;HGNC")
    expect_equal(annotated_df$EntityId[1],        "3236;3430")
    expect_equal(annotated_df$EntityName[1],      "EGFR;ERBB2")

    # A repeated member contributes its grounding once
    expect_equal(annotated_df$EntityId[2],   "3236")
    expect_equal(annotated_df$EntityName[2], "EGFR")
    expect_equal(annotated_df$IsKinase[2],   TRUE)
})

test_that("Metabolite groups keep every namespace Gilda returns per member", {
    df <- data.frame(Protein = "glucose;citrate", stringsAsFactors = FALSE)
    local_mocked_bindings(
        .callGroundEntitiesFromGildaApi = function(textInputs, keep_only = NULL, organisms = NULL) {
            list(glucose = list(ns   = c("CHEBI", "MESH"),
                                id   = c("17234", "D005947"),
                                name = c("glucose", "Glucose")),
                 citrate = list(ns = "CHEBI", id = "133748", name = "citrate"))
        }
    )
    annotated_df <- annotateProteinInfoFromIndra(df, "Metabolite")

    expect_equal(annotated_df$EntityNamespace, "CHEBI;MESH;CHEBI")
    expect_equal(annotated_df$EntityId,        "17234;D005947;133748")
    expect_equal(annotated_df$EntityName,      "glucose;Glucose;citrate")
    expect_true(is.na(annotated_df$UniprotId))
    expect_true(is.na(annotated_df$IsKinase))
})

test_that("a protein group's groundings fan out into separate query nodes", {
    df <- data.frame(Protein = "P13747;P23132", stringsAsFactors = FALSE)
    local_mocked_bindings(
        .callGetHgncIdsFromUniprotIdsApi = function(uniprotIds) {
            list(P13747 = "4931", P23132 = "10012")
        },
        .callGetHgncNamesFromHgncIdsApi = function(hgncIds) {
            list(`4931` = "HLA-E", `10012` = "RAD23A")
        },
        .callIsTranscriptionFactorApi = function(genes) list(),
        .callIsKinaseApi = function(genes) list(),
        .callIsPhosphataseApi = function(genes) list()
    )
    annotated_df <- annotateProteinInfoFromIndra(df, "Uniprot")

    pairs <- MSstatsBioNet:::.buildCogexGroundings(
        namespaces = annotated_df$EntityNamespace,
        ids        = annotated_df$EntityId,
        force_include_other = NULL
    )
    expect_setequal(
        vapply(pairs, function(p) paste(p[[1]], p[[2]], sep = ":"), character(1)),
        c("HGNC:4931", "HGNC:10012")
    )
})
