# The entity table: one row per analyte in the MSstats results. Built by
# prepare_entities(), grounded by convert_ids(), flagged by
# select_entities(), and read by get_network(). Phase 3 of the API
# refactor; internal until the end of Phase 3.

#' Input identifier systems for entities$id_type
#'
#' Which of them a backend can convert is up to its convert_ids() method.
#' @keywords internal
#' @noRd
ID_TYPES <- c("uniprot", "uniprot_mnemonic", "hgnc_symbol",
              "ensembl_protein", "ensembl_gene", "entrez",
              "chemical_name", "inchikey", "hmdb", "chebi", "chembl")

#' Fold-change columns that prepare_entities() copies to entities$logFC
#'
#' MSstats names the column after the log base it used (log2FC by default,
#' log10FC with logTrans = 10).
#' @keywords internal
#' @noRd
LOGFC_COLUMNS <- c("log2FC", "log10FC", "logFC")

#' Pattern for MSstatsPTM site identifiers
#'
#' Matches an identifier ending in one or more \code{_<residue><position>}
#' tokens, e.g. \code{P00533_S1039_S1042}. The first capture group is the
#' parent identifier, the second the site tokens.
#' @keywords internal
#' @noRd
MSSTATS_PTM_SITE_PATTERN <- "^(.*?)((?:_[A-Z][0-9]+)+)$"

#' Required entity columns and their types
#' @keywords internal
#' @noRd
REQUIRED_ENTITY_COLUMNS <- c(
    id                = "character",
    entity_type       = "character",
    id_type           = "character",
    namespace         = "character",
    entity_id         = "character",
    entity_name       = "character",
    included_in_query = "logical"
)

#' Prepare an entity table from MSstats results
#'
#' Builds the table that \code{convert_ids()} grounds,
#' \code{select_entities()} flags, and \code{get_network()} queries. It has
#' one row per analyte, and keeps every analyte, not only the significant
#' ones, so that nodes returned by a backend can be marked as measured.
#'
#' @param df output of \code{groupComparison()}'s \code{ComparisonResult}
#' table, or any table with one row per analyte.
#' @param id_column name of the column holding the analyte identifiers.
#' @param entity_type one of the entity types (\code{"protein"},
#' \code{"ptm_site"}, \code{"metabolite"}, ...), or the name of a column of
#' \code{df} holding one per row.
#' @param id_type the identifier system of \code{id_column}
#' (\code{"uniprot"}, \code{"uniprot_mnemonic"}, \code{"hgnc_symbol"},
#' \code{"chemical_name"}, ...), or the name of a column of \code{df}
#' holding one per row. For PTM sites, the identifier system of the parent
#' protein. \code{"chemical_name"} is a metabolite, lipid, or drug name,
#' common or IUPAC (e.g. \code{"glucose"}), grounded by text matching.
#' @param organism NCBI taxon ID, as a string.
#' @param label the comparison to keep, when \code{df} has a \code{Label}
#' column with more than one value.
#' @param logfc_column the fold-change column to copy to \code{logFC}.
#' \code{NULL} uses whichever of \code{log2FC}, \code{log10FC}, or
#' \code{logFC} \code{df} has. Values are copied unchanged, in the log base
#' of the input.
#' @return data.frame with columns \code{id}, \code{entity_type},
#' \code{id_type}, \code{namespace}, \code{entity_id}, \code{entity_name}
#' (all \code{NA} until \code{convert_ids()}), \code{included_in_query}
#' (\code{TRUE} until \code{select_entities()}), \code{site} and
#' \code{parent_id} (for \code{ptm_site} rows), \code{organism}, and
#' \code{logFC} and \code{adj.pvalue} when \code{df} has them.
#' @keywords internal
#' @noRd
prepare_entities <- function(df, id_column = "Protein", entity_type, id_type,
                             organism = "9606", label = NULL,
                             logfc_column = NULL) {
    df <- as.data.frame(df)
    if (!is.character(id_column) || length(id_column) != 1 ||
        !id_column %in% colnames(df)) {
        stop("id_column must name a column of df.", call. = FALSE)
    }
    if (!is.character(organism) || length(organism) != 1 ||
        !grepl("^[0-9]+$", organism)) {
        stop("organism must be an NCBI taxon ID as a string, e.g. \"9606\".",
             call. = FALSE)
    }
    df <- .select_label(df, label)
    logfc_column <- .find_logfc_column(df, logfc_column)

    ids <- as.character(df[[id_column]])
    if (anyNA(ids) || any(!nzchar(ids))) {
        stop("Column '", id_column, "' has missing or empty identifiers.",
             call. = FALSE)
    }
    duplicated_ids <- unique(ids[duplicated(ids)])
    if (length(duplicated_ids) > 0) {
        stop("Identifiers in column '", id_column, "' must be unique. ",
             "Duplicated: ", .list_values_for_message(duplicated_ids), ".",
             call. = FALSE)
    }

    entity_types <- .resolve_value_or_column(df, entity_type, ENTITY_TYPES,
                                             "entity_type")
    id_types <- .resolve_value_or_column(df, id_type, ID_TYPES, "id_type")

    n <- nrow(df)
    entities <- data.frame(
        id                = ids,
        entity_type       = entity_types,
        id_type           = id_types,
        namespace         = rep(NA_character_, n),
        entity_id         = rep(NA_character_, n),
        entity_name       = rep(NA_character_, n),
        included_in_query = rep(TRUE, n),
        site              = rep(NA_character_, n),
        parent_id         = rep(NA_character_, n),
        organism          = rep(organism, n),
        stringsAsFactors  = FALSE
    )

    is_ptm <- entities$entity_type == "ptm_site"
    if (any(is_ptm)) {
        sites <- parse_ptm_sites(entities$id[is_ptm])
        entities$site[is_ptm] <- sites$site
        entities$parent_id[is_ptm] <- sites$parent_id
        if ("GlobalProtein" %in% colnames(df)) {
            global <- as.character(df$GlobalProtein[is_ptm])
            has_global <- !is.na(global) & nzchar(global)
            entities$parent_id[is_ptm][has_global] <- global[has_global]
        }
        no_site <- is.na(entities$site[is_ptm])
        if (any(no_site)) {
            warning(sum(no_site), " ptm_site row(s) have no site in their ",
                    "identifier and are treated as their own parent: ",
                    .list_values_for_message(entities$id[is_ptm][no_site]), ".",
                    call. = FALSE)
            entities$parent_id[is_ptm][no_site] <-
                entities$id[is_ptm][no_site]
        }
    }

    if (!is.null(logfc_column)) {
        entities$logFC <- as.numeric(df[[logfc_column]])
    }
    if ("adj.pvalue" %in% colnames(df)) {
        entities$adj.pvalue <- as.numeric(df$adj.pvalue)
    }
    entities
}

#' Split PTM site identifiers into parent and site
#'
#' Each \code{";"}-separated member of an identifier (a protein group) is
#' parsed on its own. Members without a site are kept as parents.
#'
#' For example, \code{"P1_S148;P2_S5"} gives parent \code{"P1;P2"} and site
#' \code{"S148;S5"}, and \code{"P00533_S1039_S1042"} gives parent
#' \code{"P00533"} and site \code{"S1039_S1042"}.
#'
#' @param ids character vector of identifiers, e.g. \code{"P1_S148;P2_S5"}.
#' @param pattern Perl regular expression with two capture groups, the
#' parent and the site tokens. The default, \code{MSSTATS_PTM_SITE_PATTERN},
#' matches one or more trailing \code{_<residue><position>} tokens.
#' @return data.frame with columns \code{id}, \code{parent_id}, and
#' \code{site}. \code{site} is \code{"_"}-joined for several sites on one
#' protein and \code{";"}-joined across group members, and \code{NA} when no
#' member has a site. \code{parent_id} is \code{NA} when \code{site} is.
#' @keywords internal
#' @noRd
parse_ptm_sites <- function(ids, pattern = MSSTATS_PTM_SITE_PATTERN) {
    ids <- as.character(ids)
    parsed <- vapply(ids, function(id) {
        members <- strsplit(id, ";", fixed = TRUE)[[1]]
        matches <- regmatches(members,
                              regexec(pattern, members, perl = TRUE))
        has_site <- lengths(matches) == 3
        if (!any(has_site)) {
            return(c(NA_character_, NA_character_))
        }
        parents <- members
        parents[has_site] <- vapply(matches[has_site], `[`, "", 2)
        sites <- vapply(matches[has_site], function(m) sub("^_", "", m[3]), "")
        c(paste(unique(parents), collapse = ";"),
          paste(unique(sites), collapse = ";"))
    }, character(2), USE.NAMES = FALSE)
    data.frame(id = ids,
               parent_id = parsed[1, ],
               site = parsed[2, ],
               stringsAsFactors = FALSE)
}

#' Build a table of an entity table's groundings, one row per grounding
#'
#' The one place that splits the \code{";"}-joined \code{namespace},
#' \code{entity_id}, and \code{entity_name} columns. Backends read
#' groundings through it.
#'
#' @param entities entity table from \code{prepare_entities()}
#' @param namespaces namespaces to keep, e.g. \code{"HGNC"}. \code{NULL}
#' keeps all.
#' @return data.frame with one row per (entity, grounding): \code{id},
#' \code{namespace}, \code{entity_id}, \code{entity_name}. Ungrounded rows
#' are left out.
#' @keywords internal
#' @noRd
build_grounding_table <- function(entities, namespaces = NULL) {
    .validate_entities(entities)
    grounded <- !is.na(entities$namespace) & !is.na(entities$entity_id)
    entities <- entities[grounded, , drop = FALSE]
    split_namespace <- strsplit(entities$namespace, ";", fixed = TRUE)
    split_id <- strsplit(entities$entity_id, ";", fixed = TRUE)
    split_name <- strsplit(entities$entity_name, ";", fixed = TRUE)
    n_groundings <- lengths(split_namespace)
    # A name lookup that failed for every grounding leaves entity_name NA
    split_name <- lapply(seq_along(split_name), function(i) {
        if (length(split_name[[i]]) == 1 && is.na(split_name[[i]])) {
            rep(NA_character_, n_groundings[i])
        } else {
            split_name[[i]]
        }
    })
    misaligned <- n_groundings != lengths(split_id) |
        n_groundings != lengths(split_name)
    if (any(misaligned)) {
        stop("namespace, entity_id, and entity_name must have the same ",
             "number of ';'-separated values. Misaligned: ",
             .list_values_for_message(entities$id[misaligned]), ".", call. = FALSE)
    }
    long <- data.frame(
        id          = rep(entities$id, n_groundings),
        namespace   = unlist(split_namespace, use.names = FALSE),
        entity_id   = unlist(split_id, use.names = FALSE),
        entity_name = unlist(split_name, use.names = FALSE),
        stringsAsFactors = FALSE
    )
    long$entity_name[long$entity_name %in% "NA"] <- NA_character_
    if (!is.null(namespaces)) {
        long <- long[long$namespace %in% namespaces, , drop = FALSE]
    }
    rownames(long) <- NULL
    long
}

#' Flag the entities to query
#'
#' Sets \code{included_in_query} from the statistical cutoffs. Rows that
#' fail are kept, so that \code{get_network()} can still mark them as
#' measured when a backend returns them. Rows with a missing
#' \code{adj.pvalue} are not selected.
#'
#' @param entities entity table from \code{prepare_entities()}
#' @param pvalue_cutoff keep rows with \code{adj.pvalue} below this.
#' \code{NULL} applies no cutoff.
#' @param logfc_cutoff keep rows with \code{abs(logFC)} above this, on the
#' log scale of the input.
#' \code{NULL} applies no cutoff.
#' @param direction \code{"both"}, \code{"up"} (\code{logFC > 0}), or
#' \code{"down"} (\code{logFC < 0}).
#' @param include_infinite_fc whether rows with infinite \code{logFC}
#' (detected in one condition only) are selected regardless of
#' \code{pvalue_cutoff} and \code{logfc_cutoff}. \code{direction} still
#' applies.
#' @param force_include values of \code{id}, or \code{"namespace:identifier"}
#' groundings (e.g. \code{"HGNC:1234"}), selected regardless of the
#' cutoffs. Entities outside the table are added with
#' \code{get_network(include_entities = )} instead.
#' @return \code{entities} with \code{included_in_query} set, and
#' \code{user_added} (\code{TRUE} for rows selected only through
#' \code{force_include}).
#' @keywords internal
#' @noRd
select_entities <- function(entities, pvalue_cutoff = NULL,
                            logfc_cutoff = NULL,
                            direction = c("both", "up", "down"),
                            include_infinite_fc = FALSE,
                            force_include = NULL) {
    .validate_entities(entities)
    direction <- match.arg(direction)
    if (!is.null(pvalue_cutoff) &&
        (!is.numeric(pvalue_cutoff) || length(pvalue_cutoff) != 1 ||
         is.na(pvalue_cutoff))) {
        stop("pvalue_cutoff must be a single number.", call. = FALSE)
    }
    if (!is.null(logfc_cutoff) &&
        (!is.numeric(logfc_cutoff) || length(logfc_cutoff) != 1 ||
         is.na(logfc_cutoff) || logfc_cutoff < 0)) {
        stop("logfc_cutoff must be a single positive numeric value.",
             call. = FALSE)
    }
    if (!is.logical(include_infinite_fc) || length(include_infinite_fc) != 1 ||
        is.na(include_infinite_fc)) {
        stop("include_infinite_fc must be TRUE or FALSE.", call. = FALSE)
    }
    needs_logfc <- !is.null(logfc_cutoff) || direction != "both" ||
        include_infinite_fc
    .require_entity_column(entities, "adj.pvalue", !is.null(pvalue_cutoff),
                           "pvalue_cutoff")
    .require_entity_column(entities, "logFC", needs_logfc,
                           "logfc_cutoff, direction, or include_infinite_fc")

    n <- nrow(entities)
    logfc <- if ("logFC" %in% colnames(entities)) {
        entities$logFC
    } else {
        rep(NA_real_, n)
    }
    infinite <- is.infinite(logfc)
    passed <- !infinite
    if ("adj.pvalue" %in% colnames(entities)) {
        passed <- passed & !is.na(entities$adj.pvalue)
    }
    if (!is.null(pvalue_cutoff)) {
        passed <- passed & entities$adj.pvalue < pvalue_cutoff
    }
    if (!is.null(logfc_cutoff)) {
        passed <- passed & !is.na(logfc) & abs(logfc) > logfc_cutoff
    }
    if (include_infinite_fc) {
        passed <- passed | infinite
    }
    if (direction == "up") {
        passed <- passed & !is.na(logfc) & logfc > 0
    } else if (direction == "down") {
        passed <- passed & !is.na(logfc) & logfc < 0
    }

    forced <- .match_force_include(entities, force_include)
    entities$included_in_query <- passed | forced
    entities$user_added <- forced & !passed
    entities
}

#' Rows of an entity table named by force_include
#' @param entities entity table
#' @param force_include \code{id} values or \code{"namespace:identifier"}
#' strings
#' @return logical vector, one per row
#' @keywords internal
#' @noRd
.match_force_include <- function(entities, force_include) {
    if (is.null(force_include)) {
        return(rep(FALSE, nrow(entities)))
    }
    if (!is.character(force_include)) {
        stop("force_include must be a character vector.", call. = FALSE)
    }
    long <- build_grounding_table(entities)
    grounding_keys <- paste(long$namespace, long$entity_id, sep = ":")
    by_grounding <- unique(long$id[grounding_keys %in% force_include])
    matched <- entities$id %in% force_include | entities$id %in% by_grounding
    unmatched <- setdiff(force_include, c(entities$id, grounding_keys))
    if (length(unmatched) > 0) {
        message(length(unmatched), " force_include value(s) match no ",
                "entity: ", .list_values_for_message(unmatched), ". To add them as ",
                "nodes, pass them to get_network(include_entities = ).")
    }
    matched
}

#' Validate an entity table
#' @param entities entity table
#' @return \code{entities}, invisibly. Stops on the first problem.
#' @keywords internal
#' @noRd
.validate_entities <- function(entities) {
    if (!is.data.frame(entities)) {
        stop("entities must be a data.frame from prepare_entities().",
             call. = FALSE)
    }
    problems <- .check_columns(entities, REQUIRED_ENTITY_COLUMNS,
                               "entities", required = TRUE)
    if (length(problems) > 0) {
        stop(paste(problems, collapse = "; "),
             ". Build entities with prepare_entities().", call. = FALSE)
    }
    if (anyNA(entities$id) || anyDuplicated(entities$id) > 0) {
        stop("entities$id must be unique and not NA.", call. = FALSE)
    }
    if (anyNA(entities$included_in_query)) {
        stop("entities$included_in_query must not be NA.", call. = FALSE)
    }
    bad_types <- setdiff(entities$entity_type, ENTITY_TYPES)
    if (length(bad_types) > 0) {
        stop("Unknown entity_type value(s): ", .list_values_for_message(bad_types),
             ". Allowed: ", paste(ENTITY_TYPES, collapse = ", "), ".",
             call. = FALSE)
    }
    invisible(entities)
}

#' Keep one comparison of a groupComparison table
#' @param df input table
#' @param label the \code{Label} value to keep, or \code{NULL}
#' @return \code{df}, restricted to \code{label}. With no \code{label} and
#' one non-missing \code{Label} value, restricted to that value, so rows with
#' a missing \code{Label} are dropped. Unchanged when \code{df} has no
#' \code{Label} column or only missing values.
#' @keywords internal
#' @noRd
.select_label <- function(df, label) {
    has_label <- "Label" %in% colnames(df)
    if (!is.null(label)) {
        if (!is.character(label) || length(label) != 1 || is.na(label)) {
            stop("label must be a single string.", call. = FALSE)
        }
        if (!has_label) {
            stop("label was given, but df has no Label column.",
                 call. = FALSE)
        }
        if (!label %in% df$Label) {
            stop("label '", label, "' is not in df$Label. Available: ",
                 .list_values_for_message(unique(df$Label)), ".", call. = FALSE)
        }
        return(df[df$Label %in% label, , drop = FALSE])
    }
    if (has_label) {
        labels <- unique(df$Label[!is.na(df$Label)])
        if (length(labels) > 1) {
            stop("df has ", length(labels), " comparisons in its Label ",
                 "column: ", .list_values_for_message(labels), ". Choose one with ",
                 "label = .", call. = FALSE)
        }
        if (length(labels) == 1) {
            return(df[df$Label %in% labels, , drop = FALSE])
        }
    }
    df
}

#' Resolve an argument that is either a vocabulary value or a column name
#'
#' Turns \code{entity_type = "protein"} into one value per row, and
#' \code{entity_type = "kind"} into the values of column \code{kind}, checked
#' against the vocabulary.
#' @param df input table
#' @param value the argument
#' @param vocabulary allowed values
#' @param arg argument name, for messages
#' @return character vector with one value per row of \code{df}
#' @keywords internal
#' @noRd
.resolve_value_or_column <- function(df, value, vocabulary, arg) {
    if (!is.character(value) || length(value) != 1 || is.na(value)) {
        stop(arg, " must be a single string.", call. = FALSE)
    }
    if (value %in% vocabulary) {
        return(rep(value, nrow(df)))
    }
    if (!value %in% colnames(df)) {
        stop(arg, " '", value, "' is neither an allowed value nor a column ",
             "of df. Allowed values: ", paste(vocabulary, collapse = ", "),
             ".", call. = FALSE)
    }
    values <- as.character(df[[value]])
    bad <- setdiff(values, vocabulary)
    if (length(bad) > 0) {
        stop("Column '", value, "' has value(s) not allowed for ", arg, ": ",
             .list_values_for_message(bad), ". Allowed: ",
             paste(vocabulary, collapse = ", "), ".", call. = FALSE)
    }
    values
}

#' Find the fold-change column of an input table
#' @param df input table
#' @param logfc_column a column name, or \code{NULL} to look for one of
#' \code{LOGFC_COLUMNS}
#' @return the column name, or \code{NULL} when \code{df} has none
#' @keywords internal
#' @noRd
.find_logfc_column <- function(df, logfc_column) {
    if (!is.null(logfc_column)) {
        if (!is.character(logfc_column) || length(logfc_column) != 1 ||
            !logfc_column %in% colnames(df)) {
            stop("logfc_column must name a column of df.", call. = FALSE)
        }
        return(logfc_column)
    }
    found <- intersect(LOGFC_COLUMNS, colnames(df))
    if (length(found) > 1) {
        stop("df has several fold-change columns: ",
             paste(found, collapse = ", "), ". Choose one with ",
             "logfc_column = .", call. = FALSE)
    }
    if (length(found) == 0) NULL else found
}

#' Stop when a cutoff needs an entity column that is absent
#' @keywords internal
#' @noRd
.require_entity_column <- function(entities, column, needed, by) {
    if (needed && !column %in% colnames(entities)) {
        stop(by, " needs the ", column, " column in entities.", call. = FALSE)
    }
}

#' List offending values in an error, warning, or message
#'
#' Joins the values with ", " so a message can name what went wrong, e.g.
#' "Duplicated: P04637, P00533.". Shows at most \code{max_shown} values and
#' counts the rest, so a table with thousands of bad rows still gives a
#' readable message.
#' @param values the values to name
#' @param max_shown how many to show before "... (n more)"
#' @return a single string
#' @keywords internal
#' @noRd
.list_values_for_message <- function(values, max_shown = 5) {
    values <- as.character(values)
    shown <- paste(values[seq_len(min(length(values), max_shown))],
                   collapse = ", ")
    if (length(values) > max_shown) {
        shown <- paste0(shown, ", ... (", length(values) - max_shown,
                        " more)")
    }
    shown
}
