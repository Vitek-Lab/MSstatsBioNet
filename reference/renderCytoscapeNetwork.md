# Shiny render binding for cytoscapeNetwork

Shiny render binding for cytoscapeNetwork

## Usage

``` r
renderCytoscapeNetwork(expr, env = parent.frame(), quoted = FALSE)
```

## Arguments

- expr:

  An expression that generates an HTML widget (or a
  [promise](https://rstudio.github.io/promises/) of an HTML widget).

- env:

  The environment in which to evaluate `expr`.

- quoted:

  Is `expr` a quoted expression (with
  [`quote()`](https://rdrr.io/r/base/substitute.html))? This is useful
  if you want to save an expression in a variable.
