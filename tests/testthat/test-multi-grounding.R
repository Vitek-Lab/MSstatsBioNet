# Tests for the multi-grounding fan-out, membership round-trip, and the
# post-split < 400 guard introduced for the "Compound" / Entity* contract.

pair_str <- function(p) paste(p[[1]], p[[2]], sep = ":")

# ----- .buildCogexGroundings fan-out -----

test_that(".buildCogexGroundings fans out semicolon-joined (ns, id) pairs", {
    pairs <- MSstatsBioNet:::.buildCogexGroundings(
        namespaces = c("HGNC;CHEBI", "HGNC"),
        ids        = c("3815;17234", "1097"),
        force_include_other = NULL
    )
    expect_setequal(vapply(pairs, pair_str, character(1)),
                    c("HGNC:3815", "CHEBI:17234", "HGNC:1097"))
})

test_that(".buildCogexGroundings appends force_include_other groundings", {
    pairs <- MSstatsBioNet:::.buildCogexGroundings(
        namespaces = "HGNC",
        ids        = "1097",
        force_include_other = c("HGNC:9999", "CHEBI:4911")
    )
    pair_strings <- vapply(pairs, pair_str, character(1))
    expect_true("HGNC:1097"  %in% pair_strings)
    expect_true("HGNC:9999"  %in% pair_strings)
    expect_true("CHEBI:4911" %in% pair_strings)
})

test_that(".buildCogexGroundings deduplicates repeated pairs", {
    pairs <- MSstatsBioNet:::.buildCogexGroundings(
        namespaces = c("HGNC", "HGNC"),
        ids        = c("1097", "1097"),
        force_include_other = NULL
    )
    expect_equal(length(pairs), 1)
})

test_that(".buildCogexGroundings errors on mismatched per-row ns/id lengths", {
    expect_error(
        MSstatsBioNet:::.buildCogexGroundings(
            namespaces = "HGNC;CHEBI",
            ids        = "1097",
            force_include_other = NULL
        ),
        "positionally aligned"
    )
})

test_that(".buildCogexGroundings errors on bad force_include_other format", {
    expect_error(
        MSstatsBioNet:::.buildCogexGroundings(
            namespaces = "HGNC",
            ids        = "1097",
            force_include_other = "no_colon_here"
        ),
        "Invalid identifier format"
    )
})

# ----- .rowMatchesEndpoint + .addAdditionalMetadataToIndraEdge membership round-trip -----

test_that(".rowMatchesEndpoint matches a (ns, id) endpoint via membership in ;-split EntityId", {
    input <- data.frame(
        Protein         = c("FOO", "BAR"),
        EntityNamespace = c("HGNC;CHEBI", "HGNC"),
        EntityId        = c("3815;17234", "1097"),
        stringsAsFactors = FALSE
    )

    expect_equal(MSstatsBioNet:::.rowMatchesEndpoint(input, "CHEBI", "17234"),
                 c(TRUE,  FALSE))
    expect_equal(MSstatsBioNet:::.rowMatchesEndpoint(input, "HGNC",  "1097"),
                 c(FALSE, TRUE))
    # The id 17234 appears in FOO but only under namespace CHEBI, so a HGNC:17234
    # query must NOT match — namespace-awareness is the whole point.
    expect_equal(MSstatsBioNet:::.rowMatchesEndpoint(input, "HGNC",  "17234"),
                 c(FALSE, FALSE))
})

test_that(".addAdditionalMetadataToIndraEdge recovers original Protein from a multi-grounded endpoint", {
    input <- data.frame(
        Protein         = c("FOO", "BAR"),
        EntityNamespace = c("HGNC;CHEBI", "HGNC"),
        EntityId        = c("3815;17234", "1097"),
        stringsAsFactors = FALSE
    )
    edge <- list(
        source_id = "17234", source_ns = "CHEBI", source_name = "glucose",
        target_id = "1097",  target_ns = "HGNC",  target_name = "A1BG"
    )
    out <- MSstatsBioNet:::.addAdditionalMetadataToIndraEdge(edge, input)
    expect_equal(out$source_uniprot_id, "FOO")   # not "17234" or "glucose"
    expect_equal(out$target_uniprot_id, "BAR")   # not "1097" or "A1BG"
})

# ----- .constructNodesDataFrame carries entityName + entityId -----

test_that(".constructNodesDataFrame emits id, entityName, entityId, Site, logFC, adj.pvalue", {
    input <- data.frame(
        Protein         = c("FOO", "BAR"),
        EntityNamespace = c("HGNC;CHEBI", "HGNC"),
        EntityId        = c("3815;17234", "1097"),
        EntityName      = c("KIT;glucose", "A1BG"),
        Site            = c(NA_character_, NA_character_),
        log2FC          = c(1.5, -0.8),
        adj.pvalue      = c(0.01, 0.04),
        stringsAsFactors = FALSE
    )
    edges <- data.frame(source = c("FOO"), target = c("BAR"),
                        stringsAsFactors = FALSE)
    nodes <- MSstatsBioNet:::.constructNodesDataFrame(input, edges)
    expect_equal(colnames(nodes),
                 c("id", "entityName", "entityId", "Site", "logFC", "adj.pvalue"))
    expect_equal(nodes$entityName[nodes$id == "FOO"], "KIT;glucose")
    expect_equal(nodes$entityId[nodes$id == "FOO"],   "3815;17234")
    expect_equal(nodes$entityName[nodes$id == "BAR"], "A1BG")
})

# ----- < 400 guard counts post-split unique pairs -----

test_that(".validateGetSubnetworkFromIndraInput counts unique (ns, id) pairs AFTER ;-splitting", {
    # 200 rows × 2 pairs each = 400 unique pairs → fails the < 400 guard
    input_over <- data.frame(
        Protein         = paste0("P", 1:200),
        EntityNamespace = rep("HGNC;CHEBI", 200),
        EntityId        = paste0(1:200, ";C", 1:200),
        log2FC          = rep(1.0, 200),
        adj.pvalue      = rep(0.01, 200),
        stringsAsFactors = FALSE
    )
    expect_error(
        MSstatsBioNet:::.validateGetSubnetworkFromIndraInput(
            input_over, protein_level_data = NULL,
            sources_filter = NULL, force_include_other = NULL
        ),
        "less than 400 proteins"
    )

    # 200 rows × 1 pair each = 200 unique pairs → passes
    input_under <- data.frame(
        Protein         = paste0("P", 1:200),
        EntityNamespace = rep("HGNC", 200),
        EntityId        = as.character(1:200),
        log2FC          = rep(1.0, 200),
        adj.pvalue      = rep(0.01, 200),
        stringsAsFactors = FALSE
    )
    expect_silent(
        MSstatsBioNet:::.validateGetSubnetworkFromIndraInput(
            input_under, protein_level_data = NULL,
            sources_filter = NULL, force_include_other = NULL
        )
    )
})

# ----- Compound proteinIdType unit test (mocked Gilda) -----

test_that("annotateProteinInfoFromIndra with Compound mocks Gilda and skips gene-only flags", {
    df <- data.frame(Protein = c("glucose", "FOO"))
    local_mocked_bindings(
        .callGroundEntitiesFromGildaApi = function(textInputs, keep_only = NULL) {
            list(
                glucose = list(ns = "CHEBI",
                               id = "17234",
                               name = "glucose"),
                FOO     = list(ns = c("MESH", "CHEBI"),
                               id = c("3815", "17234"),
                               name = c("KIT",  "glucose"))
            )
        }
    )
    annotated_df <- annotateProteinInfoFromIndra(df, "Compound")

    expect_true(all(c("EntityNamespace", "EntityId", "EntityName") %in% colnames(annotated_df)))

    # UniprotId and gene-only flags must be NA for Compound (no API calls)
    expect_true(all(is.na(annotated_df$UniprotId)))
    expect_true(all(is.na(annotated_df$IsTranscriptionFactor)))
    expect_true(all(is.na(annotated_df$IsKinase)))
    expect_true(all(is.na(annotated_df$IsPhosphatase)))

    glucose_row <- annotated_df[annotated_df$Protein == "glucose", ]
    expect_equal(glucose_row$EntityNamespace, "CHEBI")
    expect_equal(glucose_row$EntityId,        "17234")
    expect_equal(glucose_row$EntityName,      "glucose")

    # Multi-grounded row — three Entity* columns are semicolon-joined and aligned
    foo_row <- annotated_df[annotated_df$Protein == "FOO", ]
    expect_equal(foo_row$EntityNamespace, "MESH;CHEBI")
    expect_equal(foo_row$EntityId,        "3815;17234")
    expect_equal(foo_row$EntityName,      "KIT;glucose")
})

# ----- Compound E2E test (mocked end-to-end; skipped if real fixture absent) -----

test_that("annotateProteinInfoFromIndra(Compound) -> getSubnetworkFromIndra E2E (mocked, real fixture)", {
    fixture_path <- system.file("extdata/groupComparisonModel_compound.csv",
                                package = "MSstatsBioNet")
    skip_if_not(nzchar(fixture_path) && file.exists(fixture_path),
                "Compound fixture not yet provided (see TODO-MSBio-20260528).")

    df <- data.table::fread(fixture_path)

    local_mocked_bindings(
        .callGroundEntitiesFromGildaApi = function(textInputs, keep_only = NULL) {
            result <- list()
            for (i in seq_along(textInputs)) {
                text_i <- as.character(textInputs[[i]])
                result[[text_i]] <- list(
                    ns   = "CHEBI",
                    id   = as.character(17000 + i),
                    name = paste0("compound_", i)
                )
            }
            result
        },
        .callIndraCogexApi = function(ns, ids, fio) list()
    )

    annotated <- annotateProteinInfoFromIndra(df, "Compound")
    expect_true(all(c("EntityNamespace", "EntityId", "EntityName") %in% colnames(annotated)))
    expect_true(any(grepl("CHEBI", annotated$EntityNamespace)))
    expect_true(all(is.na(annotated$UniprotId)))
    expect_true(all(is.na(annotated$IsTranscriptionFactor)))
})
