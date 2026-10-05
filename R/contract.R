#' Version of the edge and node contract checked by validate_network()
#' @keywords internal
#' @noRd
EDGE_CONTRACT_VERSION <- "1.0"

#' Normalized statement-type vocabulary for edges$interaction
#'
#' The INDRA statement types are the vocabulary. Other backends map their
#' relation labels onto it.
#' @keywords internal
#' @noRd
STATEMENT_TYPES <- c(
    # Regulation of activity or amount
    "Activation", "Inhibition", "IncreaseAmount", "DecreaseAmount",
    "Regulation", "Influence", "Gef", "Gap", "GtpActivation",
    # Modifications
    "Modification",
    "Phosphorylation", "Dephosphorylation",
    "Autophosphorylation", "Transphosphorylation",
    "Ubiquitination", "Deubiquitination",
    "Sumoylation", "Desumoylation",
    "Hydroxylation", "Dehydroxylation",
    "Acetylation", "Deacetylation",
    "Glycosylation", "Deglycosylation",
    "Farnesylation", "Defarnesylation",
    "Geranylgeranylation", "Degeranylgeranylation",
    "Palmitoylation", "Depalmitoylation",
    "Myristoylation", "Demyristoylation",
    "Ribosylation", "Deribosylation",
    "Methylation", "Demethylation",
    # Other relations
    "Complex", "Association", "Conversion", "Translocation"
)

#' Statement types whose edges are symmetric (directed = FALSE)
#' @keywords internal
#' @noRd
UNDIRECTED_STATEMENT_TYPES <- c("Complex", "Association")

#' Entity-type vocabulary for nodes$entity_type
#' @keywords internal
#' @noRd
ENTITY_TYPES <- c("protein", "gene", "transcript", "ptm_site", "metabolite",
                 "lipid", "drug", "complex", "family", "other")

#' Vocabulary for nodes$node_role
#' @keywords internal
#' @noRd
NODE_ROLES <- c("query", "forced", "mediator", "regulator", "target",
               "path_endpoint", "path_intermediate")

#' Required edge columns and their types
#' @keywords internal
#' @noRd
REQUIRED_EDGE_COLUMNS <- c(
    source         = "character",
    target         = "character",
    interaction    = "character",
    directed       = "logical",
    site           = "character",
    confidence     = "numeric",
    evidence_count = "integer",
    provenance_url = "character",
    statement_id   = "character",
    source_db      = "character",
    query_type     = "character"
)

#' Required node columns and their types
#'
#' entity_type, measured, selected and node_role become required once
#' getSubnetworkFromIndra() sets them. Until then they are checked only when
#' present (OPTIONAL_NODE_COLUMNS).
#' @keywords internal
#' @noRd
REQUIRED_NODE_COLUMNS <- c(id = "character")

#' Optional node columns, type-checked when present
#' @keywords internal
#' @noRd
OPTIONAL_NODE_COLUMNS <- c(
    entity_type        = "character",
    entity_name        = "character",
    namespace          = "character",
    entity_id          = "character",
    measured           = "logical",
    selected           = "logical",
    node_role          = "character",
    site               = "character",
    has_measured_sites = "logical",
    log2FC             = "numeric",
    adj.pvalue         = "numeric"
)

#' Validate a network against the edge and node contract
#'
#' Checks that a network returned by \code{\link{getSubnetworkFromIndra}}, or
#' built by hand from another source, has the columns, types, and
#' vocabularies that the visualization and filtering functions rely on.
#'
#' Required edge columns: \code{source}, \code{target} (both matching
#' \code{nodes$id}), \code{interaction} (an INDRA statement type such as
#' \code{"Activation"} or \code{"Complex"}), \code{directed} (logical),
#' \code{site} (PTM site on the target such as \code{"S148"}, \code{;}-joined
#' when there are several, or \code{NA}), \code{confidence} (in [0, 1], or
#' \code{NA} when the source provides no score), \code{evidence_count}
#' (whole number, at least 1), \code{provenance_url}, \code{statement_id}
#' (character), \code{source_db}, and \code{query_type}.
#'
#' Required node column: \code{id}. An \code{id} can repeat, once per PTM
#' site row of the same protein. When present, \code{entity_type},
#' \code{entity_name}, \code{namespace}, \code{entity_id}, \code{site},
#' \code{log2FC}, \code{adj.pvalue}, \code{measured}, \code{selected},
#' \code{node_role}, and \code{has_measured_sites} are type-checked.
#' Nodes with \code{measured == FALSE} must have \code{NA} statistics.
#'
#' Confidence values are comparable within one \code{source_db}, not across
#' sources.
#'
#' @param network list with \code{nodes} and \code{edges} data.frames.
#'
#' @return \code{network}, invisibly. Stops with an error listing every
#' problem found.
#'
#' @export
#'
#' @examples
#' network <- list(
#'     nodes = data.frame(id = c("CHK1_HUMAN", "CDC25A_HUMAN")),
#'     edges = data.frame(
#'         source = "CHK1_HUMAN",
#'         target = "CDC25A_HUMAN",
#'         interaction = "Phosphorylation",
#'         directed = TRUE,
#'         site = "S76",
#'         confidence = 0.99,
#'         evidence_count = 12L,
#'         provenance_url = paste0("https://db.indra.bio/statements/",
#'                                 "from_hash/-1234?format=html"),
#'         statement_id = "-1234",
#'         source_db = "INDRA",
#'         query_type = "subnetwork"
#'     )
#' )
#' validate_network(network)
#'
validate_network <- function(network) {
    if (!is.list(network) || !all(c("nodes", "edges") %in% names(network))) {
        stop("network must be a list with 'nodes' and 'edges' data.frames.")
    }
    nodes <- network$nodes
    edges <- network$edges
    if (!is.data.frame(nodes) || !is.data.frame(edges)) {
        stop("network$nodes and network$edges must be data.frames.")
    }
    problems <- c(
        .check_columns(edges, REQUIRED_EDGE_COLUMNS, "edges", required = TRUE),
        .check_columns(nodes, REQUIRED_NODE_COLUMNS, "nodes", required = TRUE),
        .check_columns(nodes, OPTIONAL_NODE_COLUMNS, "nodes", required = FALSE),
        .check_edge_values(edges),
        .check_node_values(nodes),
        .check_endpoints(nodes, edges)
    )
    if (length(problems) > 0) {
        stop("Network does not meet the edge and node contract (v",
             EDGE_CONTRACT_VERSION, "):\n",
             paste0("  - ", problems, collapse = "\n"),
             call. = FALSE)
    }
    invisible(network)
}

#' Check column presence and type
#' @param df data.frame to check
#' @param spec named character vector, column name -> type
#' @param label "edges" or "nodes", used in messages
#' @param required logical, whether missing columns are a problem
#' @return character vector of problems
#' @keywords internal
#' @noRd
.check_columns <- function(df, spec, label, required) {
    missing_cols <- setdiff(names(spec), colnames(df))
    problems <- character(0)
    if (required && length(missing_cols) > 0) {
        problems <- paste0(label, " is missing required column(s): ",
                          paste(missing_cols, collapse = ", "))
    }
    present <- intersect(names(spec), colnames(df))
    wrong_type <- present[!vapply(present, function(col) {
        .has_type(df[[col]], spec[[col]])
    }, logical(1))]
    if (length(wrong_type) > 0) {
        problems <- c(problems, paste0(
            label, "$", wrong_type, " must be ", spec[wrong_type]
        ))
    }
    problems
}

#' Check that a column has the contract type
#'
#' An all-NA logical column is accepted for any type, since that is what
#' data.frame() makes from a column of NA.
#' @keywords internal
#' @noRd
.has_type <- function(x, type) {
    if (is.logical(x) && all(is.na(x))) {
        return(TRUE)
    }
    switch(type,
           character = is.character(x),
           logical   = is.logical(x),
           numeric   = is.numeric(x),
           integer   = is.numeric(x) &&
               all(is.na(x) | x == round(x)),
           FALSE)
}

#' Check edge values against the vocabularies and ranges
#' @keywords internal
#' @noRd
.check_edge_values <- function(edges) {
    problems <- character(0)
    if (nrow(edges) == 0) {
        return(problems)
    }
    if ("interaction" %in% colnames(edges)) {
        unknown <- setdiff(unique(edges$interaction), STATEMENT_TYPES)
        if (length(unknown) > 0) {
            problems <- c(problems, paste0(
                "edges$interaction has value(s) outside the statement-type ",
                "vocabulary: ", paste(unknown, collapse = ", ")))
        }
    }
    for (col in c("source", "target", "interaction", "directed",
                  "statement_id", "provenance_url", "source_db",
                  "query_type")) {
        if (col %in% colnames(edges) && anyNA(edges[[col]])) {
            problems <- c(problems, paste0("edges$", col, " must not be NA"))
        }
    }
    if ("confidence" %in% colnames(edges) && is.numeric(edges$confidence)) {
        conf <- edges$confidence[!is.na(edges$confidence)]
        if (any(conf < 0 | conf > 1)) {
            problems <- c(problems, "edges$confidence must be in [0, 1] or NA")
        }
    }
    if ("evidence_count" %in% colnames(edges) &&
        is.numeric(edges$evidence_count) &&
        any(is.na(edges$evidence_count) | edges$evidence_count < 1)) {
        problems <- c(problems, "edges$evidence_count must be at least 1")
    }
    if ("site" %in% colnames(edges) && is.character(edges$site)) {
        sites <- edges$site[!is.na(edges$site)]
        bad <- sites[!grepl("^[A-Z][0-9]+(;[A-Z][0-9]+)*$", sites)]
        if (length(bad) > 0) {
            problems <- c(problems, paste0(
                "edges$site must look like 'S148' (';'-joined if several): ",
                paste(unique(bad), collapse = ", ")))
        }
    }
    problems
}

#' Check node values against the vocabularies and the latent-node rule
#' @keywords internal
#' @noRd
.check_node_values <- function(nodes) {
    problems <- character(0)
    if (nrow(nodes) == 0) {
        return(problems)
    }
    # nodes$id may repeat: PTM results have one row per site of a protein
    if ("id" %in% colnames(nodes) && anyNA(nodes$id)) {
        problems <- c(problems, "nodes$id must not be NA")
    }
    vocabularies <- list(entity_type = ENTITY_TYPES, node_role = NODE_ROLES)
    for (col in intersect(names(vocabularies), colnames(nodes))) {
        values <- nodes[[col]][!is.na(nodes[[col]])]
        unknown <- setdiff(unique(values), vocabularies[[col]])
        if (length(unknown) > 0) {
            problems <- c(problems, paste0(
                "nodes$", col, " has value(s) outside the vocabulary: ",
                paste(unknown, collapse = ", ")))
        }
    }
    if ("measured" %in% colnames(nodes)) {
        latent <- !is.na(nodes$measured) & !nodes$measured
        for (col in intersect(c("log2FC", "adj.pvalue"), colnames(nodes))) {
            if (any(latent & !is.na(nodes[[col]]))) {
                problems <- c(problems, paste0(
                    "nodes$", col, " must be NA for nodes with ",
                    "measured == FALSE"))
            }
        }
    }
    problems
}

#' Check that every edge endpoint is a node
#' @keywords internal
#' @noRd
.check_endpoints <- function(nodes, edges) {
    if (!all(c("source", "target") %in% colnames(edges)) ||
        !"id" %in% colnames(nodes)) {
        return(character(0))
    }
    dangling <- setdiff(unique(c(edges$source, edges$target)), nodes$id)
    if (length(dangling) == 0) {
        return(character(0))
    }
    paste0("edge endpoint(s) not found in nodes$id: ",
           paste(dangling, collapse = ", "))
}
