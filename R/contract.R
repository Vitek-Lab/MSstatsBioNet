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
INTERACTION_TYPES <- c(
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
UNDIRECTED_INTERACTION_TYPES <- c("Complex", "Association")

#' Entity-type vocabulary for nodes$entity_type
#' @keywords internal
#' @noRd
ENTITY_TYPES <- c("protein", "gene", "transcript", "ptm_site", "metabolite",
                 "lipid", "drug", "complex", "family", "other")

#' Vocabulary for nodes$node_role
#' @keywords internal
#' @noRd
NODE_ROLES <- c("passed_cutoffs", "user_added", "mediator",
               "upstream_regulator", "downstream_target",
               "path_start", "path_end", "path_intermediate")

#' Required edge columns and their types
#' @keywords internal
#' @noRd
REQUIRED_EDGE_COLUMNS <- c(
    source           = "character",
    target           = "character",
    interaction      = "character",
    directed         = "logical",
    site             = "character",
    confidence       = "numeric",
    evidence_count   = "integer",
    evidence_url     = "character",
    statement_id     = "character",
    backend_database = "character",
    query_type       = "character"
)

#' Required node columns and their types
#' @keywords internal
#' @noRd
REQUIRED_NODE_COLUMNS <- c(
    id                = "character",
    entity_type       = "character",
    entity_name       = "character",
    namespace         = "character",
    entity_id         = "character",
    measured          = "logical",
    included_in_query = "logical",
    node_role         = "character"
)

#' Optional node columns, type-checked when present
#' @keywords internal
#' @noRd
OPTIONAL_NODE_COLUMNS <- c(
    site               = "character",
    has_measured_sites = "logical",
    logFC              = "numeric",
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
#' (whole number, at least 1, not \code{NA}), \code{evidence_url},
#' \code{statement_id} (character), \code{backend_database}, and
#' \code{query_type}.
#' Edges of the symmetric statement types \code{"Complex"} and
#' \code{"Association"} must have \code{directed = FALSE}.
#'
#' Required node columns: \code{id}, \code{entity_type} (e.g.
#' \code{"protein"}, \code{"ptm_site"}, \code{"metabolite"},
#' \code{"family"}), \code{entity_name}, \code{namespace} and
#' \code{entity_id} (the grounding, \code{NA} when unknown),
#' \code{measured} (logical: the node is in the input data),
#' \code{included_in_query} (logical: the node was part of the query), and
#' \code{node_role} (why the node is in the network, e.g.
#' \code{"passed_cutoffs"} or \code{"user_added"}). An \code{id} can repeat,
#' once per PTM site row of the same protein. When present, \code{site},
#' \code{has_measured_sites}, \code{logFC}, and \code{adj.pvalue} are
#' type-checked.
#' Nodes with \code{measured == FALSE} must have \code{NA} statistics.
#'
#' Confidence values are comparable within one \code{backend_database}, not
#' across sources.
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
#'     nodes = data.frame(
#'         id = c("CHK1_HUMAN", "CDC25A_HUMAN"),
#'         entity_type = "protein",
#'         entity_name = c("CHEK1", "CDC25A"),
#'         namespace = "HGNC",
#'         entity_id = c("1925", "1725"),
#'         measured = TRUE,
#'         included_in_query = TRUE,
#'         node_role = "passed_cutoffs"
#'     ),
#'     edges = data.frame(
#'         source = "CHK1_HUMAN",
#'         target = "CDC25A_HUMAN",
#'         interaction = "Phosphorylation",
#'         directed = TRUE,
#'         site = "S76",
#'         confidence = 0.99,
#'         evidence_count = 12L,
#'         evidence_url = paste0("https://db.indra.bio/statements/",
#'                                 "from_hash/-1234?format=html"),
#'         statement_id = "-1234",
#'         backend_database = "INDRA",
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
#' @param input_dataframe data.frame to check
#' @param column_types named character vector, column name -> type
#' @param label "edges" or "nodes", used in messages
#' @param required logical, whether missing columns are a problem
#' @return character vector of problems
#' @keywords internal
#' @noRd
.check_columns <- function(input_dataframe, column_types, label, required) {
    missing_columns <- setdiff(names(column_types), colnames(input_dataframe))
    problems <- character(0)
    if (required && length(missing_columns) > 0) {
        problems <- paste0(label, " is missing required column(s): ",
                          paste(missing_columns, collapse = ", "))
    }
    present_columns <- intersect(names(column_types), colnames(input_dataframe))
    wrong_type <- present_columns[!vapply(present_columns, function(column) {
        .has_type(input_dataframe[[column]], column_types[[column]])
    }, logical(1))]
    if (length(wrong_type) > 0) {
        problems <- c(problems, paste0(
            label, "$", wrong_type, " must be ", column_types[wrong_type]
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
.has_type <- function(values, type) {
    if (is.logical(values) && all(is.na(values))) {
        return(TRUE)
    }
    switch(type,
           character = is.character(values),
           logical   = is.logical(values),
           numeric   = is.numeric(values),
           integer   = is.numeric(values) &&
               all(is.na(values) | values == round(values)),
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
        unknown <- setdiff(unique(edges$interaction), INTERACTION_TYPES)
        if (length(unknown) > 0) {
            problems <- c(problems, paste0(
                "edges$interaction has value(s) outside the statement-type ",
                "vocabulary: ", paste(unknown, collapse = ", ")))
        }
    }
    for (column in c("source", "target", "interaction", "directed",
                     "statement_id", "evidence_url", "backend_database",
                     "query_type")) {
        if (column %in% colnames(edges) && anyNA(edges[[column]])) {
            problems <- c(problems,
                          paste0("edges$", column, " must not be NA"))
        }
    }
    if (all(c("interaction", "directed") %in% colnames(edges)) &&
        is.logical(edges$directed)) {
        directed_symmetric_edges <-
            edges$interaction %in% UNDIRECTED_INTERACTION_TYPES &
            edges$directed %in% TRUE
        if (any(directed_symmetric_edges)) {
            problems <- c(problems, paste0(
                "edges$directed must be FALSE for symmetric statement types (",
                paste(UNDIRECTED_INTERACTION_TYPES, collapse = ", "), "); ",
                sum(directed_symmetric_edges),
                " row(s) have directed == TRUE"))
        }
    }
    if ("confidence" %in% colnames(edges) && is.numeric(edges$confidence)) {
        confidence <- edges$confidence[!is.na(edges$confidence)]
        if (any(confidence < 0 | confidence > 1)) {
            problems <- c(problems, "edges$confidence must be in [0, 1] or NA")
        }
    }
    if ("evidence_count" %in% colnames(edges)) {
        evidence_count <- edges$evidence_count
        # is.infinite() errors on list columns; .check_columns() reports those
        if (anyNA(evidence_count) ||
            (is.numeric(evidence_count) && any(is.infinite(evidence_count)))) {
            problems <- c(problems,
                          "edges$evidence_count must not be NA or infinite")
        }
        if (is.numeric(evidence_count) &&
            any(evidence_count[is.finite(evidence_count)] < 1)) {
            problems <- c(problems, "edges$evidence_count must be at least 1")
        }
    }
    if ("site" %in% colnames(edges) && is.character(edges$site)) {
        sites <- edges$site[!is.na(edges$site)]
        badly_formatted_sites <-
            sites[!grepl("^[A-Z][0-9]+(;[A-Z][0-9]+)*$", sites)]
        if (length(badly_formatted_sites) > 0) {
            problems <- c(problems, paste0(
                "edges$site must look like 'S148' (';'-joined if several): ",
                paste(unique(badly_formatted_sites), collapse = ", ")))
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
    for (column in intersect(names(vocabularies), colnames(nodes))) {
        values <- nodes[[column]][!is.na(nodes[[column]])]
        unknown <- setdiff(unique(values), vocabularies[[column]])
        if (length(unknown) > 0) {
            problems <- c(problems, paste0(
                "nodes$", column, " has value(s) outside the vocabulary: ",
                paste(unknown, collapse = ", ")))
        }
    }
    # A non-logical measured column is reported by .check_columns()
    if ("measured" %in% colnames(nodes) && is.logical(nodes$measured)) {
        latent <- !is.na(nodes$measured) & !nodes$measured
        statistics <- intersect(c("logFC", "adj.pvalue"), colnames(nodes))
        for (column in statistics) {
            if (any(latent & !is.na(nodes[[column]]))) {
                problems <- c(problems, paste0(
                    "nodes$", column, " must be NA for nodes with ",
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
    missing_endpoints <- setdiff(unique(c(edges$source, edges$target)),
                                 nodes$id)
    if (length(missing_endpoints) == 0) {
        return(character(0))
    }
    paste0("edge endpoint(s) not found in nodes$id: ",
           paste(missing_endpoints, collapse = ", "))
}
