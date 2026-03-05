files <- styler::style_pkg(
    transformers = styler::tidyverse_style(indent_by = 4),
    dry = "on"
)

if (any(files$changed)) {
    message("The following files need styling. Please run:")
    message("  styler::style_pkg(transformers = styler::tidyverse_style(indent_by = 4))")
    quit(status = 1)
}