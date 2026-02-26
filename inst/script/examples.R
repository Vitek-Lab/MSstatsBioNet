# =============================================================================
# cytoscapeNetwork – usage examples
# =============================================================================
# Install (once the package is built):
#   devtools::install_local("path/to/cytoscapeNetwork")
#   # or, during development:
#   devtools::load_all("path/to/cytoscapeNetwork")
# =============================================================================

library(cytoscapeNetwork)

# ── Example 2 · logFC colour gradient ───────────────────────────────────────
# Nodes coloured on a blue (down) → grey (neutral) → red (up) scale.

nodes_fc <- data.frame(
  id    = c("TP53",  "MDM2",  "CDKN1A", "BCL2",  "BAX"),
  logFC = c( 1.5,    -0.8,     2.1,     -1.9,     0.3),
  stringsAsFactors = FALSE
)

edges_fc <- data.frame(
  source      = c("TP53",   "TP53",   "MDM2",  "BCL2"),
  target      = c("MDM2",   "CDKN1A", "TP53",  "BAX"),
  interaction = c("Activation", "IncreaseAmount", "Inhibition", "Complex"),
  stringsAsFactors = FALSE
)

widget = cytoscapeNetwork(nodes_fc, edges_fc)


# ── Example 3 · PTM satellite nodes ─────────────────────────────────────────
# The `Site` column (underscore-separated) creates small circle child-nodes
# clustered around the parent protein.  Hover over edges to see overlap
# information when an edge target shares a PTM site with node data.

nodes_ptm <- data.frame(
  id    = c("EGFR",         "SRC",        "AKT1"),
  logFC = c( 1.2,            0.5,         -0.3),
  Site  = c("Y1068_Y1173",  "Y416",       NA),
  stringsAsFactors = FALSE
)

edges_ptm <- data.frame(
  source      = c("EGFR",           "SRC"),
  target      = c("SRC",            "AKT1"),
  interaction = c("Phosphorylation", "Activation"),
  site        = c("Y416",            NA),         # edge targets a specific site
  stringsAsFactors = FALSE
)

cytoscapeNetwork(nodes_ptm, edges_ptm, nodeFontSize = 14)


# ── Example 4 · HGNC labels + left-to-right layout ─────────────────────────

nodes_hgnc <- data.frame(
  id       = c("ENSG001", "ENSG002", "ENSG003"),
  hgncName = c("TP53",    "MDM2",    "CDKN1A"),
  logFC    = c( 1.0,      -0.5,       2.0),
  stringsAsFactors = FALSE
)

edges_hgnc <- data.frame(
  source      = c("ENSG001", "ENSG001"),
  target      = c("ENSG002", "ENSG003"),
  interaction = c("Activation", "IncreaseAmount"),
  stringsAsFactors = FALSE
)

cytoscapeNetwork(
  nodes_hgnc, edges_hgnc,
  displayLabelType = "hgncName",
  layoutOptions    = list(rankDir = "LR", rankSep = 120)
)


# ── Example 5 · Evidence links ───────────────────────────────────────────────
# Click an edge to open the evidence URL in a new tab.

edges_ev <- data.frame(
  source        = c("TP53", "MDM2"),
  target        = c("MDM2", "TP53"),
  interaction   = c("Activation", "Inhibition"),
  evidenceLink  = c(
    "https://www.ncbi.nlm.nih.gov/pubmed/10490031",
    "https://www.ncbi.nlm.nih.gov/pubmed/16474400"
  ),
  stringsAsFactors = FALSE
)

cytoscapeNetwork(nodes_min, edges_ev)


# ── Example 6 · Shiny integration ───────────────────────────────────────────

if (requireNamespace("shiny", quietly = TRUE)) {
  library(shiny)

  ui <- fluidPage(
    titlePanel("Protein Interaction Network"),
    sidebarLayout(
      sidebarPanel(
        sliderInput("font_size", "Node font size", min = 8, max = 24, value = 12),
        selectInput("layout_dir", "Layout direction",
                    choices = c("Top-Bottom" = "TB", "Left-Right" = "LR"),
                    selected = "TB")
      ),
      mainPanel(
        # Use the Shiny output binding
        cytoscapeNetworkOutput("network", height = "600px")
      )
    )
  )

  server <- function(input, output, session) {
    output$network <- renderCytoscapeNetwork({
      cytoscapeNetwork(
        nodes        = nodes_fc,
        edges        = edges_fc,
        nodeFontSize = input$font_size,
        layoutOptions = list(rankDir = input$layout_dir)
      )
    })
  }

  # shinyApp(ui, server)   # uncomment to launch
}


# ── Example 7 · Save to a standalone HTML file ──────────────────────────────

widget <- cytoscapeNetwork(nodes_ptm, edges_ptm)

htmlwidgets::saveWidget(
  widget,
  file       = "network.html",
  selfcontained = TRUE     # bundles all JS/CSS into one file
)
# browseURL("network.html")  # open in browser
