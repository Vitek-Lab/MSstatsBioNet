# Tests of the STRING backend (string_backend(), its convert_ids() and
# get_network() methods). STRING is never called: .post_string_request() is
# mocked with the responses saved by _fixtures/make_string_fixtures.R.

.read_string_fixture <- function(file) {
    jsonlite::fromJSON(readLines(test_path("_fixtures", "string", file)))
}

# Mocks STRING with the saved responses and records each request in the
# returned environment's `requests` (list of url and body). Waits are
# recorded in `waits` instead of slept.
.mock_string <- function(ids_response = .read_string_fixture(
                             "get_string_ids.json"),
                         env = parent.frame()) {
    calls <- new.env()
    calls$requests <- list()
    calls$waits <- numeric(0)
    local_mocked_bindings(
        .post_string_request = function(url, body) {
            calls$requests <- c(calls$requests,
                                list(list(url = url, body = body)))
            if (grepl("/version$", url)) {
                return(.read_string_fixture("version.json"))
            }
            if (grepl("/get_string_ids$", url)) {
                return(ids_response)
            }
            .read_string_fixture(paste0("network_", body$network_type,
                                        ".json"))
        },
        .wait_for_string = function(seconds) {
            calls$waits <- c(calls$waits, seconds)
        },
        .get_current_time = function() .fixed_string_time(),
        .env = env
    )
    calls
}

.fixed_string_time <- function() {
    as.POSIXct("2026-10-09 12:00:00", tz = "UTC")
}

.string_test_input <- function() {
    data.frame(Protein = c("P04637", "P04637-2", "Q00987", "P38936",
                           "P28482", "P0DTC2"),
               log2FC = c(1, 2, -1, 0.5, 3, 1),
               adj.pvalue = c(0.01, 0.02, 0.03, 0.04, 0.01, 0.5),
               stringsAsFactors = FALSE)
}

.prepare_string_entities <- function(df = .string_test_input()) {
    prepare_entities(df, entity_type = "protein", id_type = "uniprot")
}

# Entity table grounded to STRING without a request
STRING_TEST_IDS <- c(P04637 = "9606.ENSP00000269305",
                     Q00987 = "9606.ENSP00000258149",
                     P38936 = "9606.ENSP00000384849",
                     P28482 = "9606.ENSP00000215832")
STRING_TEST_NAMES <- c(P04637 = "TP53", Q00987 = "MDM2", P38936 = "CDKN1A",
                       P28482 = "MAPK1")

.grounded_string_entities <- function() {
    entities <- .prepare_string_entities()
    accession <- sub("-[0-9]+$", "", entities$id)
    grounded <- accession %in% names(STRING_TEST_IDS)
    entities$namespace[grounded] <- "STRING"
    entities$entity_id[grounded] <- STRING_TEST_IDS[accession[grounded]]
    entities$entity_name[grounded] <- STRING_TEST_NAMES[accession[grounded]]
    select_entities(entities, pvalue_cutoff = 0.05)
}

.request_urls <- function(calls) {
    vapply(calls$requests, `[[`, "", "url")
}

test_that("string_backend() defaults to physical, the current version", {
    backend <- string_backend()
    expect_s4_class(backend, "StringBackend")
    expect_s4_class(backend, "NetworkBackend")
    expect_equal(backend@network_type, "physical")
    expect_true(is.na(backend@version))
    expect_equal(backend@caller_identity, "MSstatsBioNet")
})

test_that("string_backend() rejects invalid arguments", {
    expect_error(string_backend("coexpression"), "should be one of")
    expect_error(string_backend(version = "latest"), "version")
    expect_error(string_backend(version = c("12.0", "12.5")), "version")
    expect_error(string_backend(caller_identity = ""), "caller_identity")
})

test_that("backend_capabilities() depends on the network type", {
    physical <- backend_capabilities(string_backend("physical"))
    expect_equal(physical$query_types, "subnetwork")
    expect_equal(physical$interaction_types, "Complex")
    expect_equal(physical$id_conversions,
                 list(protein = "uniprot", ptm_site = "uniprot"))
    expect_equal(physical$evidence_sources$database,
                 c("experiments", "database"))
    expect_equal(physical$evidence_sources$text_mined, "textmining")
    expect_true("coexpression" %in% physical$evidence_sources$predicted)

    regulatory <- backend_capabilities(string_backend("regulatory"))
    expect_equal(regulatory$interaction_types,
                 c("Activation", "Inhibition", "Regulation"))
    expect_length(regulatory$evidence_sources$predicted, 0)
    expect_equal(backend_capabilities(string_backend("functional"))$
                     interaction_types, "Association")
})

test_that("convert_ids() grounds accessions to STRING proteins", {
    calls <- .mock_string()
    expect_message(
        entities <- convert_ids(string_backend(), .prepare_string_entities()),
        "1 UniProt accession\\(s\\) not found in STRING: P0DTC2")
    expect_equal(entities$namespace,
                 c(rep("STRING", 5), NA))
    expect_equal(entities$entity_id[1:5],
                 unname(STRING_TEST_IDS[c("P04637", "P04637", "Q00987",
                                          "P38936", "P28482")]))
    expect_equal(entities$entity_name[1:5],
                 c("TP53", "TP53", "MDM2", "CDKN1A", "MAPK1"))
})

test_that("convert_ids() sends unique accessions without isoform suffixes", {
    calls <- .mock_string()
    suppressMessages(convert_ids(string_backend(), .prepare_string_entities()))
    first_ids_request <- calls$requests[[2]]
    expect_match(first_ids_request$url, "/get_string_ids$")
    expect_equal(strsplit(first_ids_request$body$identifiers, "\r")[[1]],
                 c("P04637", "Q00987", "P38936", "P28482", "P0DTC2"))
    expect_equal(first_ids_request$body$species, "9606")
    expect_equal(first_ids_request$body$limit, 1)
    expect_equal(first_ids_request$body$caller_identity, "MSstatsBioNet")
})

test_that("convert_ids() resends accessions STRING collapsed into one row", {
    # Two accessions of one protein: STRING answers only the first
    df <- data.frame(Protein = c("P04637", "P04637-2", "Q00987"),
                     adj.pvalue = 0.01)
    first <- data.frame(queryIndex = 0L, queryItem = "P04637",
                        stringId = "9606.ENSP00000269305",
                        ncbiTaxonId = 9606L, preferredName = "TP53")
    second <- data.frame(queryIndex = 0L, queryItem = "Q00987",
                         stringId = "9606.ENSP00000258149",
                         ncbiTaxonId = 9606L, preferredName = "MDM2")
    responses <- list(first, second)
    sent <- list()
    local_mocked_bindings(
        .post_string_request = function(url, body) {
            if (grepl("/version$", url)) {
                return(.read_string_fixture("version.json"))
            }
            sent[[length(sent) + 1]] <<- body$identifiers
            response <- responses[[1]]
            responses <<- responses[-1]
            response
        },
        .wait_for_string = function(seconds) NULL
    )
    entities <- convert_ids(string_backend(), .prepare_string_entities(df))
    expect_equal(sent, list("P04637\rQ00987", "Q00987"))
    expect_equal(entities$entity_name, c("TP53", "TP53", "MDM2"))
})

test_that("convert_ids() drops rows STRING matched to another query", {
    # A row whose queryItem differs from the accession sent at its
    # queryIndex is not trusted
    response <- data.frame(queryIndex = c(0L, 1L),
                           queryItem = c("P04637", "SOMETHING_ELSE"),
                           stringId = c("9606.ENSP00000269305", "9606.X"),
                           ncbiTaxonId = c(9606L, 9606L),
                           preferredName = c("TP53", "X"))
    .mock_string(ids_response = response)
    df <- data.frame(Protein = c("P04637", "Q00987"), adj.pvalue = 0.01)
    expect_message(
        entities <- convert_ids(string_backend(),
                                .prepare_string_entities(df)),
        "not found in STRING: Q00987")
    expect_equal(entities$entity_name, c("TP53", NA))
})

test_that("convert_ids() grounds each member of a protein group", {
    .mock_string()
    df <- data.frame(Protein = c("P04637;Q00987", "P38936"),
                     adj.pvalue = 0.01)
    entities <- convert_ids(string_backend(), .prepare_string_entities(df))
    expect_equal(entities$namespace[1], "STRING;STRING")
    expect_equal(entities$entity_name[1], "TP53;MDM2")
})

test_that("convert_ids() rejects non-accessions before any request", {
    calls <- .mock_string()
    df <- data.frame(Protein = c("P04637", "TP53_HUMAN", "TP53"),
                     adj.pvalue = 0.01)
    expect_error(convert_ids(string_backend(), .prepare_string_entities(df)),
                 "Not accessions: TP53_HUMAN, TP53")
    expect_length(calls$requests, 0)

    entities <- prepare_entities(data.frame(Protein = "TP53"),
                                 entity_type = "protein",
                                 id_type = "hgnc_symbol")
    expect_error(convert_ids(string_backend(), entities),
                 "StringBackend can't convert entity_type / id_type: protein / hgnc_symbol")
    expect_length(calls$requests, 0)
})

test_that("convert_ids() keeps INDRA's groundings and replaces STRING's", {
    .mock_string()
    entities <- .prepare_string_entities()
    entities$namespace[1] <- "HGNC"
    entities$entity_id[1] <- "11998"
    entities$entity_name[1] <- "TP53"
    entities$namespace[3] <- "HGNC;STRING"
    entities$entity_id[3] <- "6973;9606.OLD"
    entities$entity_name[3] <- "MDM2;OLD"
    entities <- suppressMessages(convert_ids(string_backend(), entities))
    expect_equal(entities$namespace[1], "HGNC;STRING")
    expect_equal(entities$entity_id[1], "11998;9606.ENSP00000269305")
    expect_equal(entities$entity_name[1], "TP53;TP53")
    expect_equal(entities$entity_id[3], "6973;9606.ENSP00000258149")
    # An unmapped row keeps nothing of STRING
    expect_true(is.na(entities$namespace[6]))
})

test_that("INDRA's convert_ids() keeps STRING's groundings", {
    local_mocked_bindings(
        .callGetHgncIdsFromUniprotIdsApi = function(uniprotIds, cogex_url) {
            list(P04637 = "11998")
        },
        .callGetHgncNamesFromHgncIdsApi = function(hgncIds, cogex_url) {
            list(`11998` = "TP53")
        }
    )
    entities <- .grounded_string_entities()
    entities <- convert_ids(indra_backend(), entities)
    expect_equal(entities$namespace[1], "STRING;HGNC")
    expect_equal(entities$entity_id[1], "9606.ENSP00000269305;11998")
    # Q00987 doesn't map to HGNC here: its STRING grounding stays
    expect_equal(entities$namespace[3], "STRING")
    expect_true(is.na(entities$namespace[6]))
})

test_that("INDRA reads only its own groundings of a shared entity table", {
    .mock_indra_response()
    indra_only <- .selected_input()
    shared <- indra_only
    shared$namespace <- ifelse(is.na(shared$namespace), "STRING",
                               paste0(shared$namespace, ";STRING"))
    shared$entity_id <- ifelse(is.na(shared$entity_id), "9606.X",
                               paste0(shared$entity_id, ";9606.X"))
    shared$entity_name <- ifelse(is.na(shared$entity_name), "X",
                                 paste0(shared$entity_name, ";X"))
    network <- suppressMessages(get_network(indra_backend(), shared))
    expected <- suppressMessages(get_network(indra_backend(), indra_only))
    expect_identical(network, expected)
    expect_false(any(grepl("STRING", network$nodes$namespace)))
})

test_that("INDRA's get_entity_properties() ignores STRING groundings", {
    local_mocked_bindings(
        .callIsKinaseApi = function(genes, cogex_url) list(MAPK1 = TRUE)
    )
    entities <- .grounded_string_entities()
    entities$namespace[5] <- "HGNC;STRING"
    entities$entity_id[5] <- "6871;9606.ENSP00000215832"
    entities$entity_name[5] <- "MAPK1;MAPK1"
    entities <- get_entity_properties(indra_backend(), entities,
                                      properties = "is_kinase")
    expect_equal(entities$is_kinase, c(NA, NA, NA, NA, TRUE, NA))
})

test_that("get_network() maps physical interactions to the contract", {
    .mock_string()
    expect_message(
        network <- get_network(string_backend(version = "12.5"),
                               .grounded_string_entities()),
        "STRING physical subnetwork: how are 5 selected proteins")
    edges <- network$edges
    expect_true(all(edges$interaction == "Complex"))
    expect_true(all(!edges$directed))
    expect_true(all(is.na(edges$site)))
    expect_true(all(is.na(edges$evidence_count)))
    expect_true(all(is.na(edges$sign)))
    expect_true(all(edges$interaction_raw == "physical"))
    expect_true(all(edges$backend_database == "STRING"))
    mdm2_tp53 <- edges[edges$source == "Q00987" & edges$target == "P04637", ]
    expect_equal(mdm2_tp53$confidence, 0.999)
    expect_equal(mdm2_tp53$evidence_sources,
                 "database;experiments;textmining")
    expect_equal(mdm2_tp53$statement_id,
                 "string:physical:9606.ENSP00000258149|9606.ENSP00000269305")
    expect_equal(mdm2_tp53$evidence_url, paste0(
        "https://version-12-5.string-db.org/cgi/network?identifiers=",
        "9606.ENSP00000258149%0d9606.ENSP00000269305&species=9606"))
    expect_no_error(validate_network(network))
})

test_that("get_network() gives isoforms sharing a STRING protein each edge", {
    .mock_string()
    network <- get_network(string_backend(version = "12.5"),
                           .grounded_string_entities())
    edges <- network$edges
    to_mdm2 <- edges[edges$source == "Q00987" &
                         edges$target %in% c("P04637", "P04637-2"), ]
    expect_equal(sort(to_mdm2$target), c("P04637", "P04637-2"))
    expect_equal(length(unique(to_mdm2$statement_id)), 1)
    tp53_nodes <- network$nodes[network$nodes$entity_name == "TP53", ]
    expect_equal(tp53_nodes$logFC, c(1, 2))
})

test_that("get_network() maps regulatory signs to interaction types", {
    .mock_string()
    network <- get_network(string_backend("regulatory", version = "12.5"),
                           .grounded_string_entities())
    edges <- network$edges
    expect_true(all(edges$directed))
    pick <- function(source, target) {
        edges[edges$source == source & edges$target == target, ]
    }
    expect_equal(pick("P28482", "P04637")$interaction, "Activation")
    expect_equal(pick("P28482", "P04637")$sign, 1L)
    expect_equal(pick("P28482", "P04637")$interaction_raw, "regulatory:pos")
    expect_equal(pick("Q00987", "P38936")$interaction, "Inhibition")
    expect_equal(pick("Q00987", "P38936")$sign, -1L)
    expect_equal(pick("Q00987", "P04637")$interaction, "Regulation")
    expect_true(is.na(pick("Q00987", "P04637")$sign))
    # Both directions of a pair are separate statements
    expect_equal(nrow(pick("P04637", "Q00987")), 1)
    expect_false(pick("P04637", "Q00987")$statement_id ==
                     pick("Q00987", "P04637")$statement_id)
    expect_equal(pick("Q00987", "P04637")$evidence_sources,
                 "database;experiments;textmining")
})

test_that("get_network() keeps an undirected pair listed in both orders once", {
    response <- .read_string_fixture("network_physical.json")
    reversed <- response[1, ]
    reversed[, c("stringId_A", "stringId_B")] <-
        reversed[, c("stringId_B", "stringId_A")]
    reversed[, c("preferredName_A", "preferredName_B")] <-
        reversed[, c("preferredName_B", "preferredName_A")]
    local_mocked_bindings(
        .post_string_request = function(url, body) rbind(response, reversed)
    )
    network <- get_network(string_backend(version = "12.5"),
                           .grounded_string_entities())
    expected <- suppressMessages({
        .mock_string()
        get_network(string_backend(version = "12.5"),
                    .grounded_string_entities())
    })
    expect_equal(nrow(network$edges), nrow(expected$edges))
})

test_that("get_network() stops when STRING returns no edges", {
    for (network_type in STRING_NETWORK_TYPES) {
        expect_equal(nrow(.parse_string_network(data.frame(), network_type)),
                     0)
    }
    local_mocked_bindings(.post_string_request = function(url, body) list())
    expect_error(get_network(string_backend("regulatory", version = "12.5"),
                             .grounded_string_entities()),
                 "No edges remain")
})

test_that("get_network() maps functional interactions to Association", {
    .mock_string()
    network <- get_network(string_backend("functional", version = "12.5"),
                           .grounded_string_entities())
    expect_true(all(network$edges$interaction == "Association"))
    expect_true(any(grepl("coexpression", network$edges$evidence_sources)))
})

test_that("get_network() sends min_confidence as required_score", {
    calls <- .mock_string()
    backend <- string_backend(version = "12.5")
    entities <- .grounded_string_entities()
    get_network(backend, entities)
    expect_equal(calls$requests[[1]]$body$required_score, 0L)
    network <- get_network(string_backend("regulatory", version = "12.5"),
                           entities, min_confidence = 0.9)
    expect_equal(calls$requests[[2]]$body$required_score, 900L)
    # Also applied after parsing (the fixture has lower scores)
    expect_true(all(network$edges$confidence >= 0.9))
})

test_that("get_network() sends the STRING proteins of the selected rows", {
    calls <- .mock_string()
    entities <- .grounded_string_entities()
    entities$included_in_query[entities$id == "P38936"] <- FALSE
    get_network(string_backend(version = "12.5"), entities)
    body <- calls$requests[[1]]$body
    expect_equal(strsplit(body$identifiers, "\r")[[1]],
                 unname(STRING_TEST_IDS[c("P04637", "Q00987", "P28482")]))
    expect_equal(body$species, "9606")
    expect_equal(body$network_type, "physical")
})

test_that("get_network() filters by evidence channel and interaction type", {
    .mock_string()
    backend <- string_backend("regulatory", version = "12.5")
    entities <- .grounded_string_entities()
    network <- get_network(backend, entities, evidence_sources = "database")
    expect_true(all(grepl("(^|;)database(;|$)", network$edges$evidence_sources)))
    expect_false(any(network$edges$source == "P38936"))

    network <- get_network(backend, entities,
                           interaction_types = "Inhibition")
    expect_equal(unique(network$edges$interaction), "Inhibition")

    expect_error(get_network(backend, entities,
                             evidence_sources = "coexpression"),
                 "must be among: experiments, database, textmining")
})

test_that("get_network() keeps NA evidence counts unless min_evidence > 1", {
    .mock_string()
    backend <- string_backend(version = "12.5")
    entities <- .grounded_string_entities()
    expect_gt(nrow(get_network(backend, entities, min_evidence = 1)$edges), 0)
    expect_message(
        expect_error(get_network(backend, entities, min_evidence = 2),
                     "No edges remain"),
        "Dropping 7 edge\\(s\\) with no evidence count")
})

test_that("get_network() adds STRING proteins from include_entities", {
    .mock_string()
    entities <- .grounded_string_entities()
    # CDKN1A becomes an entity outside the table
    entities <- entities[entities$id != "P38936", ]
    network <- get_network(string_backend(version = "12.5"), entities,
                           include_entities = "STRING:9606.ENSP00000384849")
    latent <- network$nodes[network$nodes$id == "CDKN1A", ]
    expect_equal(nrow(latent), 1)
    expect_false(latent$measured)
    expect_equal(latent$entity_type, "protein")
    expect_equal(latent$namespace, "STRING")
    expect_equal(latent$node_role, "user_added")

    expect_error(get_network(string_backend(version = "12.5"), entities,
                             include_entities = "HGNC:1234"),
                 "STRING groundings")
})

test_that("get_network() checks its input before any request", {
    calls <- .mock_string()
    backend <- string_backend(version = "12.5")
    expect_error(get_network(backend,
                             select_entities(.prepare_string_entities())),
                 "No selected entity has a STRING grounding")
    entities <- .grounded_string_entities()
    entities$organism[1] <- "10090"
    expect_error(get_network(backend, entities), "one organism")
    expect_error(get_network(backend, .grounded_string_entities(),
                             min_confidence = 2), "min_confidence")
    expect_length(calls$requests, 0)
})

test_that("get_network() records the STRING version in provenance", {
    .mock_string()
    network <- get_network(string_backend(), .grounded_string_entities())
    provenance <- network$provenance
    expect_equal(provenance$backend_database, "STRING")
    expect_equal(provenance$backend_version, "12.5")
    expect_equal(provenance$backend_url, "https://version-12-5.string-db.org")
    expect_equal(provenance$organism, "9606")
    expect_equal(provenance$retrieved_at, .fixed_string_time())
    expect_match(provenance$parameters, '"network_type":"physical"')
})

test_that("a backend resolves the current version once and keeps its host", {
    calls <- .mock_string()
    backend <- string_backend()
    entities <- suppressMessages(convert_ids(backend,
                                             .prepare_string_entities()))
    entities <- select_entities(entities, pvalue_cutoff = 0.05)
    suppressMessages(get_network(backend, entities))
    urls <- .request_urls(calls)
    expect_equal(sum(grepl("/version$", urls)), 1)
    expect_equal(urls[1], "https://string-db.org/api/json/version")
    expect_true(all(startsWith(urls[-1],
                               "https://version-12-5.string-db.org/api/json/")))
})

test_that("a pinned version needs no version request", {
    calls <- .mock_string()
    suppressMessages(get_network(string_backend(version = "12.0"),
                                 .grounded_string_entities()))
    expect_equal(.request_urls(calls),
                 "https://version-12-0.string-db.org/api/json/network")
})

test_that("STRING requests are at least one second apart", {
    clock <- as.POSIXct("2026-10-09 12:00:00", tz = "UTC")
    waits <- numeric(0)
    local_mocked_bindings(
        .post_string_request = function(url, body) list(),
        .get_string_clock = function() clock,
        .wait_for_string = function(seconds) waits <<- c(waits, seconds)
    )
    .string_request_times$last <- NULL
    .send_string_request("https://example.org/a", list())
    expect_length(waits, 0)
    clock <- clock + 0.25
    .send_string_request("https://example.org/b", list())
    expect_equal(waits, 0.75)
    clock <- clock + 5
    .send_string_request("https://example.org/c", list())
    expect_equal(waits, 0.75)
})

test_that("STRING has no evidence text or curations", {
    .mock_string()
    network <- get_network(string_backend(version = "12.5"),
                           .grounded_string_entities())
    expect_error(get_evidence(string_backend(), network$edges),
                 "StringBackend does not support get_evidence")
    expect_error(get_curations(string_backend(), network$edges),
                 "StringBackend does not support get_curations")
    expect_error(get_entity_properties(string_backend(),
                                       .grounded_string_entities()),
                 "StringBackend does not support get_entity_properties")
    # filter_by_curation() keeps edges with no evidence count, and asks no
    # backend about them
    expect_identical(filter_by_curation(network), network)
})

test_that("the STRING backend meets the shared backend contract", {
    .mock_string()
    for (network_type in STRING_NETWORK_TYPES) {
        expect_backend_contract(string_backend(network_type,
                                               version = "12.5"),
                                .grounded_string_entities())
    }
})

# An INDRA statement TP53 Complex MDM2, built from the first statement of
# the saved CoGEx response
.indra_tp53_mdm2_statement <- function() {
    statement <- readRDS(system.file("extdata/indraResponse.rds",
                                     package = "MSstatsBioNet"))[[1]]
    statement$source_ns <- "HGNC"
    statement$source_id <- "11998"
    statement$source_name <- "TP53"
    statement$target_ns <- "HGNC"
    statement$target_id <- "6973"
    statement$target_name <- "MDM2"
    statement
}

test_that("INDRA and STRING networks of one entity table share node ids", {
    .mock_string()
    local_mocked_bindings(
        .callIndraCogexApi = function(ns, ids, fio, cogex_url) {
            list(.indra_tp53_mdm2_statement())
        }
    )
    entities <- .grounded_string_entities()
    hgnc <- c(P04637 = "11998", `P04637-2` = "11998", Q00987 = "6973",
              P38936 = "1784", P28482 = "6871")
    rows <- match(names(hgnc), entities$id)
    entities$namespace[rows] <- paste0("HGNC;", entities$namespace[rows])
    entities$entity_id[rows] <- paste0(hgnc, ";", entities$entity_id[rows])
    entities$entity_name[rows] <- paste0(entities$entity_name[rows], ";",
                                         entities$entity_name[rows])

    indra <- suppressMessages(get_network(indra_backend(), entities))
    string <- suppressMessages(get_network(string_backend(version = "12.5"),
                                           entities))
    expect_setequal(c(indra$edges$source, indra$edges$target),
                    c("P04637", "P04637-2", "Q00987"))
    expect_true(all(c(indra$edges$source, indra$edges$target) %in%
                    string$nodes$id))
    expect_equal(unique(indra$nodes$namespace), "HGNC")
    expect_equal(unique(string$nodes$namespace), "STRING")

    merged <- suppressMessages(merge_networks(indra, string,
                                              entities = entities))
    expect_no_error(validate_network(merged))
    expect_setequal(merged$edges$backend_database, c("INDRA", "STRING"))
    expect_equal(nrow(merged$edges),
                 nrow(indra$edges) + nrow(string$edges))
    expect_equal(anyDuplicated(merged$nodes$id), 0)
    expect_equal(merged$provenance$backend_database, c("INDRA", "STRING"))
})

test_that("filter_by_curation() on a merged network asks only INDRA", {
    .mock_string()
    asked <- character(0)
    local_mocked_bindings(
        .callIndraCogexApi = function(ns, ids, fio, cogex_url) {
            list(.indra_tp53_mdm2_statement())
        },
        .get_incorrect_curation_count = function(statement_id, curation_url) {
            asked <<- c(asked, statement_id)
            1L
        }
    )
    entities <- .grounded_string_entities()
    entities$namespace[1] <- "STRING;HGNC"
    entities$entity_id[1] <- "9606.ENSP00000269305;11998"
    entities$entity_name[1] <- "TP53;TP53"
    entities$namespace[3] <- "STRING;HGNC"
    entities$entity_id[3] <- "9606.ENSP00000258149;6973"
    entities$entity_name[3] <- "MDM2;MDM2"
    indra <- suppressMessages(get_network(indra_backend(), entities))
    string <- suppressMessages(get_network(string_backend(version = "12.5"),
                                           entities))
    merged <- suppressMessages(merge_networks(indra, string))
    filtered <- suppressMessages(filter_by_curation(merged))
    expect_equal(unique(asked), unique(indra$edges$statement_id))
    # The INDRA edge had 1 evidence, curated as incorrect
    expect_equal(unique(filtered$edges$backend_database), "STRING")
    expect_equal(nrow(filtered$edges), nrow(string$edges))
})

test_that(".filter_by_min_evidence() keeps NA counts only at 1 or below", {
    edges <- data.frame(evidence_count = c(NA, 1L, 3L))
    expect_equal(.filter_by_min_evidence(edges, 1)$evidence_count,
                 c(NA, 1L, 3L))
    expect_equal(.filter_by_min_evidence(edges, 0)$evidence_count,
                 c(NA, 1L, 3L))
    expect_message(kept <- .filter_by_min_evidence(edges, 2),
                   "Dropping 1 edge\\(s\\) with no evidence count")
    expect_equal(kept$evidence_count, 3L)
    expect_silent(.filter_by_min_evidence(edges[2:3, , drop = FALSE], 2))
})

test_that(".replace_groundings() and .keep_groundings() round-trip", {
    entities <- .prepare_string_entities()[1:3, ]
    entities$namespace <- c("HGNC;STRING", "STRING", NA)
    entities$entity_id <- c("11998;9606.A", "9606.B", NA)
    entities$entity_name <- c("NA;TP53", "TP53", NA)
    indra_only <- .keep_groundings(entities, .is_indra_namespace)
    expect_equal(indra_only$namespace, c("HGNC", NA, NA))
    expect_true(is.na(indra_only$entity_name[1]))
    string_only <- .keep_groundings(entities, .is_string_namespace)
    expect_equal(string_only$entity_id, c("9606.A", "9606.B", NA))
    # Unchanged when every grounding is kept
    expect_identical(.keep_groundings(indra_only, .is_indra_namespace),
                     indra_only)

    replaced <- .replace_groundings(
        entities, 1:3,
        data.frame(namespace = c("STRING", NA, "STRING"),
                   entity_id = c("9606.C", NA, "9606.D"),
                   entity_name = c("C", NA, "D")),
        .is_string_namespace)
    expect_equal(replaced$namespace, c("HGNC;STRING", NA, "STRING"))
    expect_equal(replaced$entity_id, c("11998;9606.C", NA, "9606.D"))
    expect_equal(replaced$entity_name, c("NA;C", NA, "D"))
})
