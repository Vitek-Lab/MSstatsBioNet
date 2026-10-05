# ----- Entity table (Phase 3a of the API refactor) -----

.group_comparison <- function() {
    data.table::fread(
        system.file("extdata/groupComparisonModel.csv", package = "MSstatsBioNet")
    )
}

# Entities with the groundings already in groupComparisonModel.csv
.grounded_entities <- function(input = .group_comparison()) {
    entities <- prepare_entities(input, entity_type = "protein",
                                 id_type = "uniprot")
    entities$namespace <- as.character(input$EntityNamespace)
    entities$entity_id <- as.character(input$EntityId)
    entities$entity_name <- as.character(input$EntityName)
    entities
}

test_that("prepare_entities() builds one row per analyte", {
    input <- .group_comparison()
    entities <- prepare_entities(input, entity_type = "protein",
                                 id_type = "uniprot")
    expect_s3_class(entities, "data.frame")
    expect_false(data.table::is.data.table(entities))
    expect_equal(nrow(entities), nrow(input))
    expect_identical(entities$id, as.character(input$Protein))
    expect_identical(colnames(entities), c(
        "id", "entity_type", "id_type", "namespace", "entity_id",
        "entity_name", "included_in_query", "site", "parent_id", "organism",
        "log2FC", "adj.pvalue"))
    expect_true(all(entities$entity_type == "protein"))
    expect_true(all(entities$id_type == "uniprot"))
    expect_true(all(is.na(entities$namespace)))
    expect_true(all(entities$included_in_query))
    expect_true(all(is.na(entities$site)))
    expect_true(all(entities$organism == "9606"))
    expect_identical(entities$log2FC, input$log2FC)
    expect_identical(entities$adj.pvalue, input$adj.pvalue)
})

test_that("prepare_entities() leaves the caller's data.table unchanged", {
    input <- .group_comparison()
    before <- data.table::copy(input)
    prepare_entities(input, entity_type = "protein", id_type = "uniprot")
    expect_identical(input, before)
})

test_that("prepare_entities() reads id_column, organism, and omits absent statistics", {
    df <- data.frame(gene = c("TP53", "EGFR"))
    entities <- prepare_entities(df, id_column = "gene",
                                 entity_type = "protein",
                                 id_type = "hgnc_symbol", organism = "10090")
    expect_identical(entities$id, c("TP53", "EGFR"))
    expect_true(all(entities$organism == "10090"))
    expect_false(any(c("log2FC", "adj.pvalue") %in% colnames(entities)))
})

test_that("prepare_entities() takes entity_type and id_type from columns", {
    df <- data.frame(Protein = c("P04637", "glucose"),
                     kind = c("protein", "metabolite"),
                     system = c("uniprot", "name"))
    entities <- prepare_entities(df, entity_type = "kind", id_type = "system")
    expect_identical(entities$entity_type, c("protein", "metabolite"))
    expect_identical(entities$id_type, c("uniprot", "name"))
})

test_that("prepare_entities() rejects invalid arguments", {
    df <- data.frame(Protein = c("P04637", "P00533"))
    expect_error(prepare_entities(df, id_column = "Gene",
                                  entity_type = "protein", id_type = "uniprot"),
                 "id_column must name a column")
    expect_error(prepare_entities(df, entity_type = "peptide",
                                  id_type = "uniprot"),
                 "entity_type 'peptide' is neither an allowed value nor a column")
    expect_error(prepare_entities(df, entity_type = "protein",
                                  id_type = "refseq"),
                 "id_type 'refseq'")
    expect_error(prepare_entities(df, entity_type = c("protein", "gene"),
                                  id_type = "uniprot"),
                 "entity_type must be a single string")
    expect_error(prepare_entities(df, entity_type = "protein",
                                  id_type = "uniprot", organism = "human"),
                 "NCBI taxon ID")
    expect_error(prepare_entities(df, entity_type = "protein",
                                  id_type = "uniprot", organism = 9606),
                 "NCBI taxon ID")
    mixed <- data.frame(Protein = c("P04637", "x"), kind = c("protein", "rna"))
    expect_error(prepare_entities(mixed, entity_type = "kind",
                                  id_type = "uniprot"),
                 "Column 'kind' has value\\(s\\) not allowed for entity_type: rna")
})

test_that("prepare_entities() rejects missing and duplicated identifiers", {
    expect_error(prepare_entities(data.frame(Protein = c("P04637", NA)),
                                  entity_type = "protein", id_type = "uniprot"),
                 "missing or empty identifiers")
    expect_error(prepare_entities(data.frame(Protein = c("P04637", "")),
                                  entity_type = "protein", id_type = "uniprot"),
                 "missing or empty identifiers")
    expect_error(prepare_entities(data.frame(Protein = c("P04637", "P04637")),
                                  entity_type = "protein", id_type = "uniprot"),
                 "must be unique. Duplicated: P04637")
})

test_that("prepare_entities() errors on several comparisons unless label names one", {
    df <- data.frame(Protein = c("P04637", "P04637", "P00533"),
                     Label = c("A vs B", "A vs C", "A vs B"),
                     log2FC = c(1, 2, 3), adj.pvalue = c(0.01, 0.02, 0.03))
    expect_error(prepare_entities(df, entity_type = "protein",
                                  id_type = "uniprot"),
                 "df has 2 comparisons in its Label column: A vs B, A vs C")
    entities <- prepare_entities(df, entity_type = "protein",
                                 id_type = "uniprot", label = "A vs C")
    expect_identical(entities$id, "P04637")
    expect_identical(entities$log2FC, 2)
    expect_error(prepare_entities(df, entity_type = "protein",
                                  id_type = "uniprot", label = "B vs C"),
                 "label 'B vs C' is not in df\\$Label")
    expect_error(prepare_entities(df["Protein"][1, , drop = FALSE],
                                  entity_type = "protein",
                                  id_type = "uniprot", label = "A vs B"),
                 "no Label column")
    one_label <- df[df$Label == "A vs B", ]
    expect_equal(nrow(prepare_entities(one_label, entity_type = "protein",
                                       id_type = "uniprot")), 2)
})

test_that("prepare_entities() parses sites only for ptm_site rows", {
    df <- data.frame(Protein = c("P00533_S1039_S1042", "PC_A1"),
                     kind = c("ptm_site", "metabolite"))
    entities <- prepare_entities(df, entity_type = "kind", id_type = "name")
    expect_identical(entities$site, c("S1039_S1042", NA))
    expect_identical(entities$parent_id, c("P00533", NA))
    expect_identical(entities$id, c("P00533_S1039_S1042", "PC_A1"))
})

test_that("prepare_entities() takes parent_id from GlobalProtein when present", {
    df <- data.frame(Protein = c("P00533_S1064", "P04637_S15"),
                     GlobalProtein = c("P00533", NA))
    entities <- prepare_entities(df, entity_type = "ptm_site",
                                 id_type = "uniprot")
    expect_identical(entities$parent_id, c("P00533", "P04637"))
    expect_identical(entities$site, c("S1064", "S15"))
})

test_that("prepare_entities() warns on ptm_site rows without a site", {
    df <- data.frame(Protein = c("P00533_S1064", "P04637"))
    expect_warning(
        entities <- prepare_entities(df, entity_type = "ptm_site",
                                     id_type = "uniprot"),
        "1 ptm_site row\\(s\\) have no site in their identifier.*P04637")
    expect_identical(entities$parent_id, c("P00533", "P04637"))
    expect_identical(entities$site, c("S1064", NA))
})

test_that("parse_ptm_sites() splits MSstatsPTM identifiers", {
    parsed <- parse_ptm_sites(c("P00533_S1039_S1042", "P00533_S1064",
                                "CLH1_HUMAN_S148", "P1;P2_S148",
                                "P1_S148;P2_T5", "O00217", "P1_S148_extra"))
    expect_identical(colnames(parsed), c("id", "parent_id", "site"))
    expect_identical(parsed$parent_id, c("P00533", "P00533", "CLH1_HUMAN",
                                         "P1;P2", "P1;P2", NA, NA))
    expect_identical(parsed$site, c("S1039_S1042", "S1064", "S148", "S148",
                                    "S148;T5", NA, NA))
})

test_that("parse_ptm_sites() takes a custom pattern and handles empty input", {
    parsed <- parse_ptm_sites("P00533-S1064", pattern = "^(.*?)(-[A-Z][0-9]+)$")
    expect_identical(parsed$parent_id, "P00533")
    expect_identical(parsed$site, "-S1064")
    empty <- parse_ptm_sites(character(0))
    expect_equal(nrow(empty), 0)
    expect_identical(colnames(empty), c("id", "parent_id", "site"))
})

test_that("groundings() returns one row per grounding", {
    entities <- prepare_entities(data.frame(Protein = c("A", "B", "C")),
                                 entity_type = "protein", id_type = "uniprot")
    entities$namespace <- c("HGNC;CHEBI", "HGNC", NA)
    entities$entity_id <- c("1097;28748", "7715", NA)
    entities$entity_name <- c("BRAF;NA", NA, NA)
    long <- groundings(entities)
    expect_identical(long$id, c("A", "A", "B"))
    expect_identical(long$namespace, c("HGNC", "CHEBI", "HGNC"))
    expect_identical(long$entity_id, c("1097", "28748", "7715"))
    expect_identical(long$entity_name, c("BRAF", NA, NA))
    expect_identical(groundings(entities, namespaces = "CHEBI")$id, "A")
    expect_equal(nrow(groundings(entities, namespaces = "UP")), 0)
})

test_that("groundings() rejects misaligned groundings", {
    entities <- prepare_entities(data.frame(Protein = "A"),
                                 entity_type = "protein", id_type = "uniprot")
    entities$namespace <- "HGNC;HGNC"
    entities$entity_id <- "1097"
    expect_error(groundings(entities), "Misaligned: A")
})

test_that("entity validation reports missing columns and bad values", {
    entities <- prepare_entities(data.frame(Protein = c("A", "B")),
                                 entity_type = "protein", id_type = "uniprot")
    expect_error(groundings(data.frame(id = "A")),
                 "entities is missing required column\\(s\\): entity_type")
    expect_error(groundings(list(id = "A")), "must be a data.frame")
    bad <- entities
    bad$included_in_query <- c("yes", "no")
    expect_error(select_entities(bad), "included_in_query must be logical")
    bad <- entities
    bad$included_in_query[1] <- NA
    expect_error(select_entities(bad), "must not be NA")
    bad <- entities
    bad$id[2] <- "A"
    expect_error(select_entities(bad), "must be unique")
    bad <- entities
    bad$entity_type[1] <- "peptide"
    expect_error(select_entities(bad), "Unknown entity_type value\\(s\\): peptide")
})

test_that("select_entities() flags rows and drops none", {
    entities <- .grounded_entities()
    selected <- select_entities(entities, pvalue_cutoff = 0.01)
    expect_equal(nrow(selected), nrow(entities))
    expect_identical(selected$included_in_query,
                     !is.na(entities$adj.pvalue) &
                         !is.infinite(entities$log2FC) &
                         entities$adj.pvalue < 0.01)
    expect_false(any(selected$user_added))
})

test_that("select_entities() selects the rows today's input filter keeps", {
    input <- .group_comparison()
    # Add the cases the filter treats specially
    input$log2FC[1:2] <- c(Inf, -Inf)
    input$adj.pvalue[2:3] <- NA
    input$log2FC[4] <- NA
    entities <- .grounded_entities(input)
    fio <- c(paste0(input$EntityNamespace[5:6], ":", input$EntityId[5:6]),
             "HGNC:0000000")
    grid <- expand.grid(pvalue = list(NULL, 0.01), logfc = list(NULL, 1),
                        direction = c("both", "up", "down"),
                        infinite = c(FALSE, TRUE),
                        force = list(NULL, fio),
                        stringsAsFactors = FALSE)
    for (i in seq_len(nrow(grid))) {
        args <- grid[i, ]
        old <- suppressMessages(.filterGetSubnetworkFromIndraInput(
            input, args$pvalue[[1]], args$logfc[[1]], args$force[[1]],
            args$infinite, args$direction))
        new <- suppressMessages(select_entities(
            entities, pvalue_cutoff = args$pvalue[[1]],
            logfc_cutoff = args$logfc[[1]], direction = args$direction,
            include_infinite_fc = args$infinite,
            force_include = args$force[[1]]))
        expect_setequal(new$id[new$included_in_query], old$Protein)
    }
})

test_that("select_entities() marks rows selected only by force_include as user_added", {
    entities <- .grounded_entities()
    passing <- entities$id[which.min(entities$adj.pvalue)]
    failing <- entities$id[which.max(entities$adj.pvalue)]
    failing_row <- entities[entities$id == failing, ]
    grounding <- paste0(failing_row$namespace, ":", failing_row$entity_id)
    for (force in list(failing, grounding)) {
        selected <- select_entities(entities, pvalue_cutoff = 0.005,
                                    force_include = c(passing, force))
        expect_true(selected$included_in_query[selected$id == failing])
        expect_true(selected$user_added[selected$id == failing])
        expect_true(selected$included_in_query[selected$id == passing])
        expect_false(selected$user_added[selected$id == passing])
    }
})

test_that("select_entities() reports force_include values that match no entity", {
    entities <- .grounded_entities()
    expect_message(select_entities(entities, force_include = "HGNC:0000000"),
                   "1 force_include value\\(s\\) match no entity: HGNC:0000000.*include_entities")
    expect_error(select_entities(entities, force_include = 1234),
                 "force_include must be a character vector")
})

test_that("select_entities() is idempotent", {
    entities <- .grounded_entities()
    strict <- select_entities(entities, pvalue_cutoff = 0.001)
    expect_identical(select_entities(strict), select_entities(entities))
})

test_that("select_entities() works without statistics when no cutoff needs them", {
    entities <- prepare_entities(data.frame(Protein = c("A", "B")),
                                 entity_type = "protein", id_type = "uniprot")
    expect_true(all(select_entities(entities)$included_in_query))
    expect_error(select_entities(entities, pvalue_cutoff = 0.05),
                 "pvalue_cutoff needs the adj.pvalue column")
    expect_error(select_entities(entities, direction = "up"),
                 "needs the log2FC column")
})

test_that("select_entities() rejects invalid cutoffs", {
    entities <- .grounded_entities()
    expect_error(select_entities(entities, pvalue_cutoff = "0.05"),
                 "pvalue_cutoff must be a single number")
    expect_error(select_entities(entities, pvalue_cutoff = c(0.01, 0.05)),
                 "pvalue_cutoff must be a single number")
    expect_error(select_entities(entities, logfc_cutoff = -1),
                 "logfc_cutoff must be a single positive numeric value")
    expect_error(select_entities(entities, include_infinite_fc = NA),
                 "include_infinite_fc must be TRUE or FALSE")
    expect_error(select_entities(entities, direction = "sideways"),
                 "should be one of")
})
