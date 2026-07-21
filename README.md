# MSstatsBioNet

<!-- badges: start -->
[![Bioconductor Release Build](https://bioconductor.org/shields/build/release/bioc/MSstatsBioNet.svg)](https://bioconductor.org/checkResults/release/bioc-LATEST/MSstatsBioNet/)
[![Codecov test coverage](https://codecov.io/github/Vitek-Lab/MSstatsBioNet/graph/badge.svg?token=SCPSPMTOEF)](https://codecov.io/github/Vitek-Lab/MSstatsBioNet)
[![Bioconductor Downloads Rank](https://bioconductor.org/shields/downloads/release/MSstatsBioNet.svg)](https://bioconductor.org/packages/stats/bioc/MSstatsBioNet/)
[![Years in Bioconductor](https://bioconductor.org/shields/years-in-bioc/MSstatsBioNet.svg)](https://bioconductor.org/packages/release/bioc/html/MSstatsBioNet.html#since)
[![License: Artistic-2.0](https://img.shields.io/badge/license-Artistic--2.0-blue.svg)](https://opensource.org/licenses/Artistic-2.0)
<!-- badges: end -->

MSstatsBioNet is an R/Bioconductor package for network analysis and enrichment
of MSstats differential abundance results in the context of prior-knowledge
biomolecular networks. It takes the output of MSstats (or MSstatsTMT /
MSstatsPTM) differential abundance analysis, queries network databases for the
interactions among the analyzed proteins, and filters, contextualizes, and
visualizes the resulting subnetworks. Notably, it integrates with
[INDRA](https://github.com/sorgerlab/indra), a database of biological networks
assembled from the literature using text mining, enabling interpretation of
proteomic and phosphoproteomic results against past published knowledge.

MSstatsBioNet is part of the [MSstats](https://github.com/Vitek-Lab/MSstats)
family of packages, developed and maintained by the
[Vitek Lab](https://olga-vitek-lab.khoury.northeastern.edu/) at Northeastern
University. The package and its documentation are also available at
[msstats.org](http://msstats.org).

## Installation

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("MSstatsBioNet")
```

The development version can be installed directly from this repository:

```r
BiocManager::install("Vitek-Lab/MSstatsBioNet", ref = "devel")
```

## Quick Start

```r
library(MSstatsBioNet)

# Example MSstats differential abundance results (groupComparison output)
input <- data.table::fread(system.file("extdata/groupComparisonModel.csv",
                                        package = "MSstatsBioNet"))

# Retrieve the subnetwork of interactions among these proteins from INDRA
subnetwork <- getSubnetworkFromIndra(input)

head(subnetwork$nodes)
head(subnetwork$edges)

# Visualize the network (e.g. in Cytoscape, or export to HTML)
cytoscapeNetwork(subnetwork$nodes, subnetwork$edges)
```

## Input Formats

MSstatsBioNet works on **differential abundance results**, not raw search-tool
output. Its main entry point, `getSubnetworkFromIndra()`, accepts the
`ComparisonResult` table produced by the group-comparison functions across the
MSstats ecosystem:

| Upstream package | Function producing input |
| --- | --- |
| [MSstats](https://github.com/Vitek-Lab/MSstats) | `groupComparison()` |
| [MSstatsTMT](https://github.com/Vitek-Lab/MSstatsTMT) | `groupComparisonTMT()` |
| [MSstatsPTM](https://github.com/Vitek-Lab/MSstatsPTM) | `groupComparisonPTM()` |

The input table provides, per protein and comparison, the log2 fold change,
p-value, and adjusted p-value used for filtering and network coloring. UniProt
identifiers can be annotated with `annotateProteinInfoFromIndra()`.

- **Databases supported:** INDRA
- **Filtering options:** p-value filter, context/topic-based filtering
  (`filterSubnetworkByContext()`)
- **Visualization options:** Cytoscape Desktop (`cytoscapeNetwork()`), in-browser
  preview (`previewNetworkInBrowser()`), standalone HTML export
  (`exportNetworkToHTML()`), and Shiny integration
  (`cytoscapeNetworkOutput()` / `renderCytoscapeNetwork()`)

## Documentation

- [MSstatsBioNet overview](vignettes/MSstatsBioNet.Rmd) — getting started
- [Cytoscape visualization](vignettes/Cytoscape-Visualization.Rmd)
- [Filter by context](vignettes/Filter-By-Context.Rmd)
- [PTM analysis](vignettes/PTM-Analysis.Rmd)
- [Official website: msstats.org](http://msstats.org)
- [Bioconductor package page and reference manual](https://bioconductor.org/packages/MSstatsBioNet)

## Getting Help / Reporting Bugs

- **Questions about usage, statistical methods, or troubleshooting:** please
  post to the [MSstats Google Group](https://groups.google.com/forum/#!forum/msstats).
  This is monitored by the development team and searchable, so it's the fastest
  way to get help and to see if your question has already been answered.
- **Bug reports and feature requests for this repository:** please open a
  [GitHub issue](https://github.com/Vitek-Lab/MSstatsBioNet/issues).

## References

If you use MSstatsBioNet, please cite:

1. Wu A, Kohler D, Navada P, Robbins J, Boyle G, Boshart A, Karis K, Neefjes J,
   Konvalinka A, Sarthy J, Pino L, Gyori B, Vitek O. **MSstatsBioNet: Integrating
   Statistical Analyses with Prior Knowledge Biomolecular Networks for
   Quantitative Proteomics and Phosphoproteomics.** *bioRxiv*. 2026.
   [DOI: 10.64898/2026.07.09.737605](https://doi.org/10.64898/2026.07.09.737605)

## Funding

MSstats development has been supported by the Chan Zuckerberg Initiative's
[Essential Open Source Software for Science](https://chanzuckerberg.com/eoss/proposals/).

## License

MSstatsBioNet is released under the [Artistic-2.0](https://opensource.org/licenses/Artistic-2.0)
license. However, its dependencies may have different licenses. Notably, INDRA
is distributed under the [BSD 2-Clause](https://opensource.org/license/bsd-2-clause)
license, and INDRA's knowledge sources may have different licenses for commercial
applications. Please refer to the
[INDRA README](https://github.com/sorgerlab/indra?tab=readme-ov-file#indra-modules)
for more information on its knowledge sources and their associated licenses.
