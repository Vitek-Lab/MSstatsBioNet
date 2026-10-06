# ----- INDRA backend and get_network() (Phase 2 of the API refactor) -----

.selected_input <- function() {
    input <- data.table::fread(
        system.file("extdata/groupComparisonModel.csv", package = "MSstatsBioNet")
    )
    .filterGetSubnetworkFromIndraInput(input, pvalueCutoff = NULL,
                                       logfc_cutoff = NULL,
                                       force_include_other = NULL,
                                       include_infinite_fc = FALSE,
                                       direction = "both")
}

.mock_indra_response <- function(env = parent.frame()) {
    local_mocked_bindings(
        .callIndraCogexApi = function(ns, ids, fio, cogex_url) {
            readRDS(system.file("extdata/indraResponse.rds", package = "MSstatsBioNet"))
        },
        .env = env
    )
}

test_that("indra_backend() defaults to the public CoGEx and Gilda URLs", {
    backend <- indra_backend()
    expect_s4_class(backend, "IndraBackend")
    expect_s4_class(backend, "NetworkBackend")
    expect_equal(backend@cogex_url, "https://discovery.indra.bio")
    expect_equal(backend@grounding_url, "https://grounding.indra.bio")
})

test_that("indra_backend() rejects an invalid cogex_url", {
    expect_error(indra_backend(cogex_url = character(0)), "cogex_url")
    expect_error(indra_backend(cogex_url = c("a", "b")), "cogex_url")
    expect_error(indra_backend(cogex_url = ""), "cogex_url")
    expect_error(indra_backend(cogex_url = NA_character_), "cogex_url")
})

test_that("indra_backend() rejects an invalid grounding_url", {
    expect_error(indra_backend(grounding_url = character(0)), "grounding_url")
    expect_error(indra_backend(grounding_url = ""), "grounding_url")
    expect_error(indra_backend(grounding_url = NA_character_), "grounding_url")
})

test_that("get_network() reproduces the pinned golden output", {
    .mock_indra_response()
    network <- get_network(indra_backend(), .selected_input(), subnetwork_query())
    golden <- readRDS(test_path("_fixtures", "golden_subnetwork.rds"))
    expect_identical(network, golden)
})

test_that("get_network() runs subnetwork_query() when query is missing", {
    .mock_indra_response()
    input <- .selected_input()
    expect_identical(
        get_network(indra_backend(), input),
        get_network(indra_backend(), input, subnetwork_query())
    )
})

test_that("get_network() passes the arguments on to the old filters", {
    .mock_indra_response()
    input <- .selected_input()
    network <- get_network(indra_backend(), input,
                           interaction_types = "Complex", min_evidence = 2)
    expect_true(all(network$edges$interaction == "Complex"))
    expect_true(all(network$edges$evidence_count >= 2))

    network <- get_network(indra_backend(), input,
                           interaction_types = c("Activation", "Phosphorylation"))
    expect_equal(nrow(network$nodes), 3)
    expect_equal(nrow(network$edges), 2)

    network <- get_network(indra_backend(), input, evidence_sources = "signor")
    expect_gt(nrow(network$edges), 0)
    expect_true(all(grepl('"signor"', network$edges$evidence_sources)))
})

test_that("get_network() sends the backend's cogex_url", {
    sent_url <- NULL
    local_mocked_bindings(
        .callIndraCogexApi = function(ns, ids, fio, cogex_url) {
            sent_url <<- cogex_url
            readRDS(system.file("extdata/indraResponse.rds", package = "MSstatsBioNet"))
        }
    )
    get_network(indra_backend(cogex_url = "https://cogex.example.org"),
                .selected_input())
    expect_equal(sent_url, "https://cogex.example.org")
})

test_that(".callIndraCogexApi() posts to the subnetwork endpoint under cogex_url", {
    posted_url <- NULL
    local_mocked_bindings(
        POST = function(url, ...) {
            posted_url <<- url
            structure(list(), class = "response")
        },
        content = function(x, ...) list()
    )
    .callIndraCogexApi("HGNC", "1925", NULL, "https://cogex.example.org")
    expect_equal(posted_url,
                 "https://cogex.example.org/api/indra_subnetwork_relations")
})

test_that("get_network() errors for an unsupported backend and query pair", {
    where <- environment()
    setClass("UnsupportedQuery", contains = "NetworkQuery", where = where)
    on.exit(removeClass("UnsupportedQuery", where = where), add = TRUE)
    expect_error(
        get_network(indra_backend(), .selected_input(),
                    new("UnsupportedQuery")),
        "IndraBackend does not support UnsupportedQuery"
    )
})

test_that("get_network() validates its input before calling INDRA", {
    local_mocked_bindings(
        .callIndraCogexApi = function(ns, ids, fio, cogex_url) {
            stop("INDRA should not be called")
        }
    )
    input <- .selected_input()
    expect_error(get_network(indra_backend(), input[0, ]),
                 "at least one protein")
    no_entity_id <- as.data.frame(input)
    no_entity_id$EntityId <- NULL
    expect_error(get_network(indra_backend(), no_entity_id),
                 "missing required column\\(s\\): EntityId")
    expect_error(get_network(indra_backend(), input, evidence_sources = 1),
                 "evidence_sources must be a character vector")
    # The old function's sources_filter is passed on as evidence_sources
    expect_error(getSubnetworkFromIndra(input, sources_filter = 1),
                 "evidence_sources must be a character vector")
})

# ----- convert_ids() and get_annotations() (Phase 3b of the API refactor) -----

.entities <- function(ids, entity_type = "protein", id_type = "uniprot", ...) {
    suppressWarnings(prepare_entities(
        data.frame(Protein = ids, stringsAsFactors = FALSE),
        entity_type = entity_type, id_type = id_type, ...))
}

# Records every call in `calls` so tests can check inputs and URLs
.mock_cogex_id_mapping <- function(calls, env = parent.frame()) {
    local_mocked_bindings(
        .callGetUniprotIdsFromUniprotMnemonicIdsApi = function(uniprotMnemonicIds, cogex_url) {
            calls$mnemonic <- c(calls$mnemonic, list(list(ids = unlist(uniprotMnemonicIds), url = cogex_url)))
            list(CLH1_HUMAN = "Q00610", P53_HUMAN = "P04637")
        },
        .callGetHgncIdsFromUniprotIdsApi = function(uniprotIds, cogex_url) {
            calls$hgnc <- c(calls$hgnc, list(list(ids = unlist(uniprotIds), url = cogex_url)))
            list(Q00610 = "2092", P04637 = "11998", P13747 = "4962", P23132 = "10012")
        },
        .callGetHgncNamesFromHgncIdsApi = function(hgncIds, cogex_url) {
            calls$names <- c(calls$names, list(list(ids = unlist(hgncIds), url = cogex_url)))
            list(`2092` = "CLTC", `11998` = "TP53", `4962` = "HLA-E", `10012` = "RAD23A")
        },
        .env = env
    )
}

.mock_gilda <- function(calls, env = parent.frame()) {
    local_mocked_bindings(
        .callGroundEntitiesFromGildaApi = function(textInputs, keep_only = NULL,
                                                   organisms = NULL, grounding_url) {
            calls$gilda <- c(calls$gilda, list(list(texts = unlist(textInputs),
                keep_only = keep_only, organisms = organisms, url = grounding_url)))
            list(TP53 = list(ns = "HGNC", id = "11998", name = "TP53"),
                 glucose = list(ns = c("CHEBI", "CHEBI"),
                                id = c("CHEBI:17234", "CHEBI:4167"),
                                name = c("glucose", "D-glucopyranose")))
        },
        .env = env
    )
}

test_that("convert_ids() grounds UniProt IDs to HGNC, pooling protein groups", {
    calls <- new.env()
    .mock_cogex_id_mapping(calls)
    entities <- convert_ids(indra_backend(cogex_url = "https://cogex.example.org"),
                            .entities(c("Q00610", "P13747;P23132", "XXXX99")))

    expect_equal(entities$namespace, c("HGNC", "HGNC;HGNC", NA))
    expect_equal(entities$entity_id, c("2092", "4962;10012", NA))
    expect_equal(entities$entity_name, c("CLTC", "HLA-E;RAD23A", NA))
    expect_setequal(calls$hgnc[[1]]$ids, c("Q00610", "P13747", "P23132", "XXXX99"))
    expect_equal(calls$hgnc[[1]]$url, "https://cogex.example.org")
    expect_equal(calls$names[[1]]$url, "https://cogex.example.org")
    expect_length(calls$hgnc, 1)
})

test_that("convert_ids() grounds PTM sites by their parent protein", {
    calls <- new.env()
    .mock_cogex_id_mapping(calls)
    entities <- convert_ids(indra_backend(),
                            .entities(c("P04637_S15", "P04637_S20", "Q00610"),
                                      entity_type = "ptm_site"))

    expect_setequal(calls$hgnc[[1]]$ids, c("P04637", "Q00610"))
    expect_equal(entities$entity_id, c("11998", "11998", "2092"))
    expect_equal(entities$site, c("S15", "S20", NA))
})

test_that("convert_ids() maps UniProt mnemonics to UniProt IDs first", {
    calls <- new.env()
    .mock_cogex_id_mapping(calls)
    entities <- convert_ids(indra_backend(),
                            .entities(c("CLH1_HUMAN", "NOPE_HUMAN"),
                                      id_type = "uniprot_mnemonic"))

    expect_equal(calls$mnemonic[[1]]$ids, c("CLH1_HUMAN", "NOPE_HUMAN"))
    expect_equal(calls$hgnc[[1]]$ids, "Q00610")
    expect_equal(entities$entity_name, c("CLTC", NA))
})

test_that("convert_ids() grounds HGNC symbols through Gilda with the entities' organism", {
    calls <- new.env()
    .mock_gilda(calls)
    entities <- convert_ids(indra_backend(grounding_url = "https://gilda.example.org"),
                            .entities(c("TP53", "NOTAGENE"), id_type = "hgnc_symbol"))

    expect_equal(entities$namespace, c("HGNC", NA))
    expect_equal(entities$entity_id, c("11998", NA))
    expect_equal(calls$gilda[[1]]$keep_only, "HGNC")
    expect_equal(calls$gilda[[1]]$organisms, list("9606"))
    expect_equal(calls$gilda[[1]]$url, "https://gilda.example.org")
})

test_that("convert_ids() grounds chemical names through Gilda without restrictions", {
    calls <- new.env()
    .mock_gilda(calls)
    entities <- convert_ids(indra_backend(),
                            .entities("glucose", entity_type = "metabolite",
                                      id_type = "chemical_name"))

    expect_equal(entities$namespace, "CHEBI;CHEBI")
    expect_equal(entities$entity_id, "CHEBI:17234;CHEBI:4167")
    expect_equal(entities$entity_name, "glucose;D-glucopyranose")
    expect_null(calls$gilda[[1]]$keep_only)
    expect_null(calls$gilda[[1]]$organisms)
})

test_that("convert_ids() makes one batch of calls per id_type in a mixed table", {
    calls <- new.env()
    .mock_cogex_id_mapping(calls)
    .mock_gilda(calls)
    input <- data.frame(Protein = c("Q00610", "glucose", "P04637", "TP53"),
                        kind = c("protein", "metabolite", "protein", "protein"),
                        system = c("uniprot", "chemical_name", "uniprot", "hgnc_symbol"))
    entities <- convert_ids(indra_backend(),
                            prepare_entities(input, entity_type = "kind",
                                             id_type = "system"))

    expect_length(calls$hgnc, 1)
    expect_length(calls$gilda, 2)
    expect_equal(entities$entity_name,
                 c("CLTC", "glucose;D-glucopyranose", "TP53", "TP53"))
    expect_equal(entities$id, input$Protein)
})

test_that("convert_ids() makes no calls for an empty table", {
    local_mocked_bindings(
        .callGetHgncIdsFromUniprotIdsApi = function(uniprotIds, cogex_url) stop("called")
    )
    entities <- convert_ids(indra_backend(), .entities(character(0)))
    expect_equal(nrow(entities), 0)
})

test_that("convert_ids() rejects id types and organisms INDRA can't convert", {
    expect_error(
        convert_ids(indra_backend(),
                    .entities("glucose", entity_type = "metabolite")),
        "can't convert entity_type / id_type: metabolite / uniprot"
    )
    expect_error(
        convert_ids(indra_backend(),
                    .entities("ENSG1", entity_type = "gene", id_type = "ensembl_gene")),
        "gene / ensembl_gene"
    )
    expect_error(
        convert_ids(indra_backend(), .entities("Q00610", organism = "10090")),
        "human .*Got organism: 10090"
    )
    expect_error(convert_ids(indra_backend(), data.frame(id = "a")),
                 "prepare_entities")
})

test_that("convert_ids() accepts a non-human organism for chemicals", {
    calls <- new.env()
    .mock_gilda(calls)
    entities <- convert_ids(indra_backend(),
                            .entities("glucose", entity_type = "metabolite",
                                      id_type = "chemical_name", organism = "10090"))
    expect_equal(entities$namespace, "CHEBI;CHEBI")
})

.annotation_entities <- function() {
    entities <- .entities(c("P04637", "P00533_S1039", "P1;P2", "glucose", "XXXX99"),
                          entity_type = "protein")
    entities$entity_type[2] <- "ptm_site"
    entities$entity_type[4] <- "metabolite"
    entities$namespace <- c("HGNC", "HGNC", "HGNC;HGNC", "CHEBI", NA)
    entities$entity_id <- c("11998", "3236", "1;2", "CHEBI:17234", NA)
    entities$entity_name <- c("TP53", "EGFR", "A;B", "glucose", NA)
    entities
}

.mock_annotation_apis <- function(calls, env = parent.frame()) {
    record <- function(field, value) {
        force(field)
        function(genes, cogex_url) {
            calls[[field]] <- list(genes = unlist(genes), url = cogex_url)
            stats::setNames(as.list(rep(value, length(genes))), unlist(genes))
        }
    }
    local_mocked_bindings(
        .callIsTranscriptionFactorApi = record("tf", TRUE),
        .callIsKinaseApi = record("kinase", FALSE),
        .callIsPhosphataseApi = record("phosphatase", FALSE),
        .env = env
    )
}

test_that("get_annotations() annotates rows with a single HGNC grounding", {
    calls <- new.env()
    .mock_annotation_apis(calls)
    entities <- get_annotations(indra_backend(cogex_url = "https://cogex.example.org"),
                                .annotation_entities())

    expect_equal(entities$is_transcription_factor, c(TRUE, TRUE, NA, NA, NA))
    expect_equal(entities$is_kinase, c(FALSE, FALSE, NA, NA, NA))
    expect_equal(entities$is_phosphatase, c(FALSE, FALSE, NA, NA, NA))
    expect_equal(calls$tf$genes, c("TP53", "EGFR"))
    expect_equal(calls$kinase$url, "https://cogex.example.org")
})

test_that("get_annotations() adds only the requested fields", {
    calls <- new.env()
    .mock_annotation_apis(calls)
    entities <- get_annotations(indra_backend(), .annotation_entities(),
                                fields = "is_kinase")

    expect_true("is_kinase" %in% colnames(entities))
    expect_false("is_transcription_factor" %in% colnames(entities))
    expect_null(calls$tf)
})

test_that("get_annotations() leaves NA where the API gives no answer", {
    local_mocked_bindings(
        .callIsTranscriptionFactorApi = function(genes, cogex_url) NULL,
        .callIsKinaseApi = function(genes, cogex_url) list(TP53 = NULL, EGFR = TRUE),
        .callIsPhosphataseApi = function(genes, cogex_url) list()
    )
    entities <- get_annotations(indra_backend(), .annotation_entities())

    expect_equal(entities$is_transcription_factor, rep(NA, 5))
    expect_equal(entities$is_kinase, c(NA, TRUE, NA, NA, NA))
    expect_equal(entities$is_phosphatase, rep(NA, 5))
})

test_that("get_annotations() makes no calls when no row can be annotated", {
    local_mocked_bindings(
        .callIsKinaseApi = function(genes, cogex_url) stop("called")
    )
    entities <- get_annotations(indra_backend(), .annotation_entities()[4:5, ],
                                fields = "is_kinase")
    expect_equal(entities$is_kinase, c(NA, NA))
})

test_that("get_annotations() rejects unknown fields", {
    expect_error(get_annotations(indra_backend(), .annotation_entities(),
                                 fields = c("is_kinase", "is_gpcr")),
                 "IndraBackend does not support annotation field\\(s\\): is_gpcr")
    expect_error(get_annotations(indra_backend(), .annotation_entities(),
                                 fields = 1),
                 "fields must be a character vector")
})

test_that("convert_ids() and get_annotations() error for a backend without them", {
    where <- environment()
    setClass("BareBackend", contains = "NetworkBackend", where = where)
    on.exit(removeClass("BareBackend", where = where), add = TRUE)
    expect_error(convert_ids(new("BareBackend"), .entities("Q00610")),
                 "BareBackend does not support convert_ids\\(\\)")
    expect_error(get_annotations(new("BareBackend"), .entities("Q00610")),
                 "BareBackend does not support get_annotations\\(\\)")
})

test_that("annotateProteinInfoFromIndra() grounds through convert_ids() and get_annotations()", {
    local_mocked_bindings(
        convert_ids = function(backend, entities, ...) {
            expect_s4_class(backend, "IndraBackend")
            expect_equal(entities$id, c("P04637", "Q00610"))
            expect_equal(entities$id_type, c("uniprot", "uniprot"))
            entities$namespace <- "HGNC"
            entities$entity_id <- c("11998", "2092")
            entities$entity_name <- c("TP53", "CLTC")
            entities
        },
        get_annotations = function(backend, entities, ...) {
            entities$is_transcription_factor <- c(TRUE, FALSE)
            entities$is_kinase <- FALSE
            entities$is_phosphatase <- FALSE
            entities
        }
    )
    df <- data.frame(Protein = c("P04637_S15", "Q00610", "P04637_S20", NA))
    annotated <- annotateProteinInfoFromIndra(df, "Uniprot")

    expect_equal(annotated$EntityName, c("TP53", "CLTC", "TP53", NA))
    expect_equal(annotated$IsTranscriptionFactor, c(TRUE, FALSE, TRUE, NA))
    expect_equal(colnames(annotated),
                 c("Protein", "GlobalProtein", "UniprotId", "EntityNamespace",
                   "EntityId", "EntityName", "IsTranscriptionFactor",
                   "IsKinase", "IsPhosphatase"))
})
