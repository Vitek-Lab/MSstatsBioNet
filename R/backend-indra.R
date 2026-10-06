#' Base URL of INDRA CoGEx
#' @keywords internal
#' @noRd
INDRA_API_URL <- "https://discovery.indra.bio"

#' Base URL of Gilda, INDRA's grounding service
#' @keywords internal
#' @noRd
GILDA_API_URL <- "https://grounding.indra.bio"

#' Create an INDRA backend
#'
#' Internal until the entity model is added (Phase 3 of the API refactor).
#' @param cogex_url base URL of INDRA CoGEx
#' @param grounding_url base URL of Gilda
#' @return an \code{IndraBackend} object
#' @importFrom methods new
#' @keywords internal
#' @noRd
indra_backend <- function(cogex_url = INDRA_API_URL,
                          grounding_url = GILDA_API_URL) {
    new("IndraBackend", cogex_url = cogex_url, grounding_url = grounding_url)
}

#' INDRA subnetwork query
#'
#' Sends the selected rows' groundings to CoGEx
#' \code{indra_subnetwork_relations}, filters the statements, and normalizes
#' them to the edge contract. Nodes are the selected rows that an edge
#' reaches, plus any \code{include_entities} endpoints.
#' @keywords internal
#' @noRd
setMethod("get_network", signature("IndraBackend", "SubnetworkQuery"),
    function(backend, entities, query, interaction_types = NULL,
             min_evidence = 1, evidence_sources = NULL,
             include_entities = NULL, ...) {
        .validateIndraSubnetworkInput(entities, evidence_sources,
                                      include_entities)
        res <- .callIndraCogexApi(entities$EntityNamespace, entities$EntityId,
                                  include_entities, backend@cogex_url)
        res <- .filterIndraResponse(res, interaction_types, min_evidence,
                                    evidence_sources)
        edges <- .constructEdgesDataFrame(res, entities)
        edges <- .filterEdgesDataFrame(edges)
        network <- list(nodes = .constructNodesDataFrame(entities, edges),
                        edges = edges)
        validate_network(network)
        network
    })

#' Validate the input of the INDRA subnetwork query
#' @param input annotated groupComparison table of the selected rows
#' @param evidence_sources evidence sources filter
#' @param force_include_other character vector of identifiers to include in
#' the network
#' @keywords internal
#' @noRd
.validateIndraSubnetworkInput <- function(input, evidence_sources, force_include_other) {
    required_cols <- c("Protein", "log2FC", "adj.pvalue",
                       "EntityNamespace", "EntityId", "EntityName")
    missing_cols <- setdiff(required_cols, colnames(input))
    if (length(missing_cols) > 0) {
        stop("Invalid Input Error: input is missing required column(s): ",
             paste(missing_cols, collapse = ", "), ".")
    }
    ids_split <- unlist(strsplit(as.character(input$EntityId),        ";"), use.names = FALSE)
    nss_split <- unlist(strsplit(as.character(input$EntityNamespace), ";"), use.names = FALSE)
    unique_pairs <- unique(paste(nss_split, ids_split, sep = ":"))
    num_proteins = length(unique_pairs) +
        ifelse(!is.null(force_include_other), length(force_include_other), 0)
    if (num_proteins >= 400) {
        stop("Invalid Input Error: INDRA query must contain less than 400 proteins.  Consider lowering your p-value cutoff")
    }
    if (nrow(input) == 0) {
        stop("Invalid Input Error: Input must contain at least one protein after filtering.")
    }
    if (!is.null(evidence_sources)) {
        if (!is.character(evidence_sources)) {
            stop("evidence_sources must be a character vector")
        }
    }
}

#' Split a protein group into its member identifiers
#'
#' A \code{Protein} value may name a protein group -- several identifiers
#' for the same quantified analyte joined by \code{";"}. Splits on
#' \code{";"}, trims surrounding whitespace and drops empty members, so a
#' plain single identifier comes back as a length-one vector.
#'
#' @param x A length-one character value, possibly \code{NA}.
#' @return A character vector of member identifiers, empty when the input
#'         holds none.
#' @keywords internal
#' @noRd
.splitProteinGroup <- function(x) {
        if (length(x) == 0 || is.na(x)) {
                return(character(0))
        }
        members <- trimws(unlist(strsplit(as.character(x), ";", fixed = TRUE),
                                 use.names = FALSE))
        return(members[nzchar(members)])
}

#' Join protein group members back into a single value
#'
#' @param members A character vector of member identifiers.
#' @return The members joined by \code{";"}, or \code{NA} when empty.
#' @keywords internal
#' @noRd
.joinProteinGroup <- function(members) {
        if (length(members) == 0) {
                return(NA_character_)
        }
        return(paste(members, collapse = ";"))
}

#' Identifier systems the INDRA backend converts, by entity type
#'
#' Proteins (and the parent proteins of PTM sites) are grounded to HGNC
#' through CoGEx or Gilda. Chemical names are grounded through Gilda to
#' whichever namespace it returns (CHEBI, PUBCHEM, CHEMBL, ...).
#' @keywords internal
#' @noRd
INDRA_ID_CONVERSIONS <- list(
    protein    = c("uniprot", "uniprot_mnemonic", "hgnc_symbol"),
    ptm_site   = c("uniprot", "uniprot_mnemonic", "hgnc_symbol"),
    metabolite = "chemical_name",
    lipid      = "chemical_name",
    drug       = "chemical_name"
)

#' Entity properties of the INDRA backend
#'
#' \code{api} names the CoGEx call, which takes a list of HGNC gene symbols.
#' A PTM site gets the properties of its parent protein.
#' @keywords internal
#' @noRd
INDRA_ENTITY_PROPERTIES <- list(
    is_transcription_factor = list(api = ".callIsTranscriptionFactorApi",
                                   entity_types = c("protein", "ptm_site")),
    is_kinase               = list(api = ".callIsKinaseApi",
                                   entity_types = c("protein", "ptm_site")),
    is_phosphatase          = list(api = ".callIsPhosphataseApi",
                                   entity_types = c("protein", "ptm_site"))
)

#' INDRA identifier conversion
#'
#' Groups rows by \code{id_type} and makes one batch of calls per group:
#' \code{uniprot} through CoGEx's UniProt-to-HGNC mapping,
#' \code{uniprot_mnemonic} through CoGEx's mnemonic-to-UniProt mapping
#' first, \code{hgnc_symbol} through Gilda restricted to HGNC and the
#' row's organism, and \code{chemical_name} through Gilda with no namespace
#' restriction. Each \code{";"}-separated member of an identifier (a protein
#' group) is grounded on its own, and the groundings are pooled onto the
#' row.
#' @keywords internal
#' @noRd
setMethod("convert_ids", "IndraBackend",
    function(backend, entities, ...) {
        .validate_entities(entities)
        .check_indra_id_conversions(entities)
        members <- .get_grounding_inputs(entities)
        for (id_type in unique(entities$id_type)) {
            rows <- which(entities$id_type == id_type)
            groundings <- switch(id_type,
                uniprot = .ground_uniprot_with_cogex(
                    members[rows], backend@cogex_url),
                uniprot_mnemonic = .ground_uniprot_with_cogex(
                    .map_uniprot_mnemonics(members[rows], backend@cogex_url),
                    backend@cogex_url),
                hgnc_symbol = .ground_text_with_gilda(
                    members[rows], backend@grounding_url, keep_only = "HGNC",
                    organisms = as.list(unique(.get_entity_organisms(entities)[rows]))),
                chemical_name = .ground_text_with_gilda(
                    members[rows], backend@grounding_url))
            entities$namespace[rows] <- groundings$namespace
            entities$entity_id[rows] <- groundings$entity_id
            entities$entity_name[rows] <- groundings$entity_name
        }
        entities
    })

#' INDRA entity properties
#'
#' Supports the properties in \code{INDRA_ENTITY_PROPERTIES}. They are
#' looked up by gene symbol, so only rows with a single HGNC grounding get
#' values; rows with several groundings (a protein group, or an ambiguous
#' name) are \code{NA}.
#' @keywords internal
#' @noRd
setMethod("get_entity_properties", "IndraBackend",
    function(backend, entities, properties = NULL, ...) {
        .validate_entities(entities)
        properties <- .resolve_entity_properties(
            backend, properties, names(INDRA_ENTITY_PROPERTIES))
        # A ";"-joined namespace (several groundings) never equals "HGNC"
        single_hgnc <- entities$namespace %in% "HGNC" &
            !is.na(entities$entity_name)
        for (property in properties) {
            spec <- INDRA_ENTITY_PROPERTIES[[property]]
            queried <- single_hgnc & entities$entity_type %in% spec$entity_types
            values <- rep(NA, nrow(entities))
            genes <- unique(entities$entity_name[queried])
            if (length(genes) > 0) {
                # Looked up by name at call time, so tests can mock the call
                call_api <- get(spec$api, mode = "function")
                response <- call_api(as.list(genes), backend@cogex_url)
                for (gene in names(response)) {
                    if (!is.null(response[[gene]])) {
                        values[queried & entities$entity_name == gene] <-
                            response[[gene]]
                    }
                }
            }
            entities[[property]] <- values
        }
        entities
    })

#' Stop when the INDRA backend can't convert some rows
#'
#' Checks each row's (\code{entity_type}, \code{id_type}) pair against
#' \code{INDRA_ID_CONVERSIONS}, and that proteins are human: CoGEx maps
#' UniProt IDs to HGNC, and Gilda is restricted to HGNC, which covers human
#' genes only.
#' @param entities entity table
#' @keywords internal
#' @noRd
.check_indra_id_conversions <- function(entities) {
    pairs <- unique(entities[, c("entity_type", "id_type")])
    supported <- vapply(seq_len(nrow(pairs)), function(i) {
        pairs$id_type[i] %in% INDRA_ID_CONVERSIONS[[pairs$entity_type[i]]]
    }, logical(1))
    if (any(!supported)) {
        unsupported <- paste(pairs$entity_type[!supported],
                             pairs$id_type[!supported], sep = " / ")
        allowed <- unlist(lapply(names(INDRA_ID_CONVERSIONS), function(type) {
            paste(type, INDRA_ID_CONVERSIONS[[type]], sep = " / ")
        }))
        stop("IndraBackend can't convert entity_type / id_type: ",
             .list_values_for_message(unsupported), ". Supported: ",
             paste(allowed, collapse = ", "), ".", call. = FALSE)
    }
    is_protein <- entities$entity_type %in% c("protein", "ptm_site")
    organisms <- unique(.get_entity_organisms(entities)[is_protein])
    non_human <- setdiff(organisms, "9606")
    if (length(non_human) > 0) {
        stop("IndraBackend grounds proteins to HGNC, which covers human ",
             "(organism \"9606\") only. Got organism: ",
             .list_values_for_message(non_human), ".", call. = FALSE)
    }
}

#' The organism of each entity row
#' @param entities entity table
#' @return character vector, \code{"9606"} where the table has no
#' \code{organism} column
#' @keywords internal
#' @noRd
.get_entity_organisms <- function(entities) {
    if ("organism" %in% colnames(entities)) {
        as.character(entities$organism)
    } else {
        rep("9606", nrow(entities))
    }
}

#' The identifiers to ground for each entity row
#' @param entities entity table
#' @return list with one character vector per row: the members of
#' \code{id}, or of \code{parent_id} for \code{ptm_site} rows
#' @keywords internal
#' @noRd
.get_grounding_inputs <- function(entities) {
    inputs <- entities$id
    if ("parent_id" %in% colnames(entities)) {
        use_parent <- entities$entity_type == "ptm_site" &
            !is.na(entities$parent_id)
        inputs[use_parent] <- entities$parent_id[use_parent]
    }
    lapply(inputs, .splitProteinGroup)
}

#' Map UniProt mnemonics to UniProt IDs through CoGEx
#' @param members list of character vectors of mnemonics, one per row
#' @param cogex_url base URL of INDRA CoGEx
#' @return list of character vectors of UniProt IDs, one per row. Members
#' that don't map are dropped.
#' @keywords internal
#' @noRd
.map_uniprot_mnemonics <- function(members, cogex_url) {
    mnemonics <- unique(unlist(members, use.names = FALSE))
    if (length(mnemonics) == 0) {
        return(members)
    }
    mapping <- .callGetUniprotIdsFromUniprotMnemonicIdsApi(as.list(mnemonics),
                                                           cogex_url)
    lapply(members, function(row_members) {
        unique(as.character(unlist(mapping[row_members], use.names = FALSE)))
    })
}

#' Ground UniProt IDs to HGNC through CoGEx
#' @param members list of character vectors of UniProt IDs, one per row
#' @param cogex_url base URL of INDRA CoGEx
#' @return data.frame with one row per element of \code{members}:
#' \code{namespace}, \code{entity_id}, \code{entity_name}, \code{";"}-joined
#' and \code{NA} for rows that don't ground. \code{entity_name} has
#' \code{"NA"} where one of several name lookups failed, and is \code{NA}
#' when all of them did.
#' @keywords internal
#' @noRd
.ground_uniprot_with_cogex <- function(members, cogex_url) {
    groundings <- .build_empty_groundings(length(members))
    uniprot_ids <- unique(unlist(members, use.names = FALSE))
    if (length(uniprot_ids) == 0) {
        return(groundings)
    }
    hgnc_mapping <- .callGetHgncIdsFromUniprotIdsApi(as.list(uniprot_ids),
                                                     cogex_url)
    hgnc_ids <- unique(as.character(unlist(hgnc_mapping, use.names = FALSE)))
    name_mapping <- list()
    if (length(hgnc_ids) > 0) {
        name_response <- .callGetHgncNamesFromHgncIdsApi(as.list(hgnc_ids),
                                                         cogex_url)
        if (!is.null(name_response)) {
            name_mapping <- name_response
        }
    }
    for (i in seq_along(members)) {
        entity_ids <- unique(as.character(
            unlist(hgnc_mapping[members[[i]]], use.names = FALSE)))
        if (length(entity_ids) == 0) {
            next
        }
        entity_names <- vapply(name_mapping[entity_ids], function(name) {
            if (is.null(name)) NA_character_ else as.character(name)[1]
        }, character(1), USE.NAMES = FALSE)
        groundings$namespace[i] <-
            .joinProteinGroup(rep("HGNC", length(entity_ids)))
        groundings$entity_id[i] <- .joinProteinGroup(entity_ids)
        if (!all(is.na(entity_names))) {
            groundings$entity_name[i] <- .joinProteinGroup(entity_names)
        }
    }
    groundings
}

#' Ground names through Gilda
#' @param members list of character vectors of names, one per row
#' @param grounding_url base URL of Gilda
#' @param keep_only namespace to keep, or \code{NULL} to keep all
#' @param organisms list of NCBI taxon IDs to restrict Gilda to, or
#' \code{NULL}
#' @return data.frame like \code{.ground_uniprot_with_cogex()}'s. Every
#' candidate Gilda returns is kept, in its ranking order, without duplicate
#' (namespace, identifier) pairs.
#' @keywords internal
#' @noRd
.ground_text_with_gilda <- function(members, grounding_url, keep_only = NULL,
                                    organisms = NULL) {
    groundings <- .build_empty_groundings(length(members))
    texts <- unique(unlist(members, use.names = FALSE))
    if (length(texts) == 0) {
        return(groundings)
    }
    grounding_map <- .callGroundEntitiesFromGildaApi(
        as.list(texts), keep_only = keep_only, organisms = organisms,
        grounding_url = grounding_url)
    if (is.null(grounding_map)) {
        return(groundings)
    }
    for (i in seq_along(members)) {
        namespaces <- character(0)
        entity_ids <- character(0)
        entity_names <- character(0)
        for (candidates in grounding_map[members[[i]]]) {
            if (is.null(candidates)) {
                next
            }
            stopifnot(length(candidates$ns) == length(candidates$id),
                      length(candidates$ns) == length(candidates$name))
            namespaces <- c(namespaces, as.character(candidates$ns))
            entity_ids <- c(entity_ids, as.character(candidates$id))
            entity_names <- c(entity_names, as.character(candidates$name))
        }
        keep <- !duplicated(paste(namespaces, entity_ids, sep = ":"))
        if (!any(keep)) {
            next
        }
        groundings$namespace[i] <- .joinProteinGroup(namespaces[keep])
        groundings$entity_id[i] <- .joinProteinGroup(entity_ids[keep])
        groundings$entity_name[i] <- .joinProteinGroup(entity_names[keep])
    }
    groundings
}

#' Ungrounded grounding columns for n rows
#' @keywords internal
#' @noRd
.build_empty_groundings <- function(n) {
    data.frame(namespace = rep(NA_character_, n),
               entity_id = rep(NA_character_, n),
               entity_name = rep(NA_character_, n),
               stringsAsFactors = FALSE)
}

#' Check the properties asked of get_entity_properties()
#' @param backend the backend, for the message
#' @param properties requested properties, or \code{NULL} for all
#' @param supported the properties the backend supports
#' @return the properties to add
#' @keywords internal
#' @noRd
.resolve_entity_properties <- function(backend, properties, supported) {
    if (is.null(properties)) {
        return(supported)
    }
    if (!is.character(properties) || anyNA(properties)) {
        stop("properties must be a character vector.", call. = FALSE)
    }
    unknown <- setdiff(properties, supported)
    if (length(unknown) > 0) {
        stop(class(backend), " does not support these entity properties: ",
             .list_values_for_message(unknown), ". Supported: ",
             paste(supported, collapse = ", "), ".", call. = FALSE)
    }
    unique(properties)
}
