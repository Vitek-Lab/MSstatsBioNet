# Tests for the multi-grounding fan-out, membership round-trip, and the
# post-split < 400 guard introduced for the "Metabolite" / Entity* contract.

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

# ----- Grounding lookup + .addAdditionalMetadataToIndraEdge membership round-trip -----

.multi_grounded_entities <- function(included_in_query = c(TRUE, TRUE)) {
    entities <- MSstatsBioNet:::.build_entities_from_annotated_input(data.frame(
        Protein         = c("FOO", "BAR"),
        log2FC          = c(1.5, -0.8),
        adj.pvalue      = c(0.01, 0.04),
        EntityNamespace = c("HGNC;CHEBI", "HGNC"),
        EntityId        = c("3815;17234", "1097"),
        EntityName      = c("KIT;glucose", "A1BG"),
        stringsAsFactors = FALSE
    ))
    entities$included_in_query <- included_in_query
    entities
}

test_that(".find_entity_rows_for_grounding matches a (namespace, identifier) grounding via membership in ;-split entity_id", {
    grounding_lookup <- MSstatsBioNet:::.build_grounding_lookup(.multi_grounded_entities())

    expect_equal(MSstatsBioNet:::.find_entity_rows_for_grounding(grounding_lookup, "CHEBI", "17234"), 1L)
    expect_equal(MSstatsBioNet:::.find_entity_rows_for_grounding(grounding_lookup, "HGNC", "1097"), 2L)
    # The id 17234 appears in FOO but only under namespace CHEBI, so a HGNC:17234
    # query must NOT match — namespace-awareness is the whole point.
    expect_equal(MSstatsBioNet:::.find_entity_rows_for_grounding(grounding_lookup, "HGNC", "17234"),
                 integer(0))
})

test_that(".find_entity_rows_for_grounding prefers rows in the query", {
    entities <- .multi_grounded_entities(c(TRUE, FALSE))
    entities$entity_id[2] <- "3815"
    grounding_lookup <- MSstatsBioNet:::.build_grounding_lookup(entities)
    expect_equal(MSstatsBioNet:::.find_entity_rows_for_grounding(grounding_lookup, "HGNC", "3815"), 1L)
    # A grounding no queried row has is matched against all rows
    entities$entity_id[2] <- "1097"
    grounding_lookup <- MSstatsBioNet:::.build_grounding_lookup(entities)
    expect_equal(MSstatsBioNet:::.find_entity_rows_for_grounding(grounding_lookup, "HGNC", "1097"), 2L)
})

test_that(".addAdditionalMetadataToIndraEdge recovers original Protein from a multi-grounded source", {
    grounding_lookup <- MSstatsBioNet:::.build_grounding_lookup(.multi_grounded_entities())
    edge <- list(
        source_id = "17234", source_ns = "CHEBI", source_name = "glucose",
        target_id = "1097",  target_ns = "HGNC",  target_name = "A1BG"
    )
    edge_with_node_ids <- MSstatsBioNet:::.addAdditionalMetadataToIndraEdge(edge, grounding_lookup)
    expect_equal(edge_with_node_ids$source_node_id, "FOO") # not "17234" or "glucose"
    expect_equal(edge_with_node_ids$target_node_id, "BAR") # not "1097" or "A1BG"
})

# ----- .build_network_nodes carries entity_name + entity_id -----

test_that(".build_network_nodes emits the node contract columns", {
    grounding_lookup <- MSstatsBioNet:::.build_grounding_lookup(.multi_grounded_entities())
    edges <- data.frame(source = c("FOO"), target = c("BAR"),
                        stringsAsFactors = FALSE)
    nodes <- MSstatsBioNet:::.build_network_nodes(
        grounding_lookup, edges, MSstatsBioNet:::.list_backend_nodes(list(), grounding_lookup))
    expect_equal(colnames(nodes),
                 c("id", "entity_type", "entity_name", "namespace", "entity_id",
                   "measured", "included_in_query", "node_role", "site",
                   "has_measured_sites", "logFC", "adj.pvalue"))
    expect_equal(nodes$entity_name[nodes$id == "FOO"], "KIT;glucose")
    expect_equal(nodes$entity_id[nodes$id == "FOO"],   "3815;17234")
    expect_equal(nodes$entity_name[nodes$id == "BAR"], "A1BG")
    expect_equal(nodes$entity_type, c("protein", "protein"))
})

# ----- < 400 guard counts post-split unique pairs -----

test_that(".validateIndraSubnetworkInput counts unique (namespace, identifier) groundings after ;-splitting", {
    .build_query_groundings <- function(namespace, entity_id) {
        entities <- MSstatsBioNet:::.build_entities_from_annotated_input(data.frame(
            Protein         = paste0("P", 1:200),
            log2FC          = rep(1.0, 200),
            adj.pvalue      = rep(0.01, 200),
            EntityNamespace = namespace,
            EntityId        = entity_id,
            EntityName      = gsub("[0-9]+", "name", entity_id),
            stringsAsFactors = FALSE
        ))
        MSstatsBioNet:::.get_groundings_to_query(entities)
    }
    # 200 rows × 2 pairs each = 400 unique pairs → fails the < 400 guard
    expect_error(
        MSstatsBioNet:::.validateIndraSubnetworkInput(
            .build_query_groundings(rep("HGNC;CHEBI", 200), paste0(1:200, ";C", 1:200)),
            evidence_sources = NULL, include_entities = NULL
        ),
        "less than 400 proteins"
    )

    # 200 rows × 1 pair each = 200 unique pairs → passes
    expect_silent(
        MSstatsBioNet:::.validateIndraSubnetworkInput(
            .build_query_groundings(rep("HGNC", 200), as.character(1:200)),
            evidence_sources = NULL, include_entities = NULL
        )
    )
})

# ----- Metabolite proteinIdType unit test (mocked Gilda) -----

test_that("annotateProteinInfoFromIndra with Metabolite mocks Gilda and skips gene-only flags", {
    df <- data.frame(Protein = c("glucose", "FOO"))
    local_mocked_bindings(
        .callGroundEntitiesFromGildaApi = function(textInputs, keep_only = NULL, organisms = NULL, grounding_url) {
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
    annotated_df <- annotateProteinInfoFromIndra(df, "Metabolite")

    expect_true(all(c("EntityNamespace", "EntityId", "EntityName") %in% colnames(annotated_df)))

    # UniprotId and gene-only flags must be NA for Metabolite (no API calls)
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

# ----- Metabolite E2E test (mocked end-to-end; skipped if real fixture absent) -----

test_that("annotateProteinInfoFromIndra(Metabolite) -> getSubnetworkFromIndra E2E (mocked, real fixture)", {
    fixture_path <- system.file("extdata/groupComparisonModel_compound.csv",
                                package = "MSstatsBioNet")
    skip_if_not(nzchar(fixture_path) && file.exists(fixture_path),
                "Metabolite fixture not yet provided (see TODO-MSBio-20260528).")

    df <- data.table::fread(fixture_path)

    local_mocked_bindings(
        .callGroundEntitiesFromGildaApi = function(textInputs, keep_only = NULL, organisms = NULL, grounding_url) {
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
        .callIndraCogexApi = function(ns, ids, fio, cogex_url) list()
    )

    annotated <- annotateProteinInfoFromIndra(df, "Metabolite")
    expect_true(all(c("EntityNamespace", "EntityId", "EntityName") %in% colnames(annotated)))
    expect_true(any(grepl("CHEBI", annotated$EntityNamespace)))
    expect_true(all(is.na(annotated$UniprotId)))
    expect_true(all(is.na(annotated$IsTranscriptionFactor)))
})
