library(mockery)

make_nodes <- function() {
    data.frame(
        id    = c("P53_HUMAN", "MDM2_HUMAN"),
        logFC = c(1.5, -1.0),
        stringsAsFactors = FALSE
    )
}

make_edges <- function() {
    data.frame(
        source      = "P53_HUMAN",
        target      = "MDM2_HUMAN",
        interaction = "Activation",
        stringsAsFactors = FALSE
    )
}

test_that("exportNetworkToHTML calls saveWidget with the correct filename", {
    tmp <- tempfile(fileext = ".html")
    
    save_widget_mock <- mock()
    stub(exportNetworkToHTML, "htmlwidgets::saveWidget", save_widget_mock)
    
    exportNetworkToHTML(make_nodes(), make_edges(), filename = tmp)
    
    expect_called(save_widget_mock, 1)
    call_args <- mock_args(save_widget_mock)[[1]]
    expect_equal(call_args$file, tmp)
})

test_that("exportNetworkToHTML calls saveWidget with selfcontained = TRUE", {
    save_widget_mock <- mock()
    stub(exportNetworkToHTML, "htmlwidgets::saveWidget", save_widget_mock)
    
    exportNetworkToHTML(make_nodes(), make_edges(), filename = tempfile(fileext = ".html"))
    
    call_args <- mock_args(save_widget_mock)[[1]]
    expect_true(call_args$selfcontained)
})

test_that("exportNetworkToHTML passes a cytoscapeNetwork widget to saveWidget", {
    save_widget_mock <- mock()
    stub(exportNetworkToHTML, "htmlwidgets::saveWidget", save_widget_mock)
    
    exportNetworkToHTML(make_nodes(), make_edges(), filename = tempfile(fileext = ".html"))
    
    widget_arg <- mock_args(save_widget_mock)[[1]][[1]]
    expect_s3_class(widget_arg, "htmlwidget")
    expect_s3_class(widget_arg, "cytoscapeNetwork")
})

test_that("exportNetworkToHTML passes nodeFontSize through to the widget", {
    save_widget_mock <- mock()
    stub(exportNetworkToHTML, "htmlwidgets::saveWidget", save_widget_mock)
    
    exportNetworkToHTML(make_nodes(), make_edges(),
                        filename      = tempfile(fileext = ".html"),
                        nodeFontSize  = 18)
    
    widget_arg <- mock_args(save_widget_mock)[[1]][[1]]
    expect_equal(widget_arg$x$node_font_size, 18)
})

test_that("exportNetworkToHTML passes displayLabelType through to the widget", {
    save_widget_mock <- mock()
    stub(exportNetworkToHTML, "htmlwidgets::saveWidget", save_widget_mock)
    
    nodes_hgnc <- make_nodes()
    nodes_hgnc$hgncName <- c("TP53", "MDM2")
    
    exportNetworkToHTML(nodes_hgnc, make_edges(),
                        filename         = tempfile(fileext = ".html"),
                        displayLabelType = "hgncName")
    
    widget_arg   <- mock_args(save_widget_mock)[[1]][[1]]
    protein_nodes <- Filter(function(el) !is.null(el$data$node_type) &&
                                el$data$node_type == "protein",
                            widget_arg$x$elements)
    labels <- sapply(protein_nodes, function(el) el$data$label)
    expect_length(labels, 2)
    expect_setequal(labels, c("TP53", "MDM2"))
})

# =============================================================================
# previewNetworkInBrowser — mock saveWidget and browseURL
# =============================================================================

test_that("previewNetworkInBrowser calls exportNetworkToHTML with a .html temp path", {
    export_mock <- mock(invisible(NULL))
    stub(previewNetworkInBrowser, "exportNetworkToHTML", export_mock)
    stub(previewNetworkInBrowser, "interactive",         mock(FALSE))
    
    previewNetworkInBrowser(make_nodes(), make_edges())
    
    expect_called(export_mock, 1)
    call_args <- mock_args(export_mock)[[1]]
    expect_true(grepl("\\.html$", call_args$filename))
})

test_that("previewNetworkInBrowser calls browseURL when interactive() is TRUE", {
    browse_mock  <- mock(invisible(NULL))
    export_mock  <- mock(invisible(NULL))
    mock_interactive <- mock(TRUE)
    
    stub(previewNetworkInBrowser, "exportNetworkToHTML", export_mock)
    stub(previewNetworkInBrowser, "browseURL",    browse_mock)
    stub(previewNetworkInBrowser, "interactive",         mock_interactive)
    
    previewNetworkInBrowser(make_nodes(), make_edges())
    
    expect_called(browse_mock, 1)
})

test_that("previewNetworkInBrowser passes the temp filename to browseURL", {
    captured_filename <- NULL
    browse_mock <- mock(invisible(NULL))
    
    stub(previewNetworkInBrowser, "exportNetworkToHTML",
         function(nodes, edges, filename, ...) {
             captured_filename <<- filename
             invisible(filename)
         })
    stub(previewNetworkInBrowser, "browseURL", browse_mock)
    stub(previewNetworkInBrowser, "interactive",      mock(TRUE))
    
    result <- previewNetworkInBrowser(make_nodes(), make_edges())
    
    browse_arg <- mock_args(browse_mock)[[1]][[1]]
    expect_equal(browse_arg, captured_filename)
    expect_equal(result,     captured_filename)
})

test_that("previewNetworkInBrowser does NOT call browseURL when non-interactive", {
    browse_mock <- mock(invisible(NULL))
    export_mock <- mock(invisible(NULL))
    
    
    stub(previewNetworkInBrowser, "exportNetworkToHTML", export_mock)
    stub(previewNetworkInBrowser, "utils::browseURL",    browse_mock)
    stub(previewNetworkInBrowser, "interactive",         mock(FALSE))
    
    previewNetworkInBrowser(make_nodes(), make_edges())
    
    expect_called(browse_mock, 0)
})

test_that("previewNetworkInBrowser passes nodeFontSize and displayLabelType through", {
    export_mock <- mock(invisible(NULL))
    stub(previewNetworkInBrowser, "exportNetworkToHTML", export_mock)
    stub(previewNetworkInBrowser, "interactive",         mock(FALSE))
    
    nodes_hgnc           <- make_nodes()
    nodes_hgnc$hgncName  <- c("TP53", "MDM2")
    
    previewNetworkInBrowser(nodes_hgnc, make_edges(),
                            displayLabelType = "hgncName",
                            nodeFontSize     = 16)
    
    call_args <- mock_args(export_mock)[[1]]
    expect_equal(call_args$displayLabelType, "hgncName")
    expect_equal(call_args$nodeFontSize,     16)
})