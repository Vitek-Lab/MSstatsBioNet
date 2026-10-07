# MSstatsBioNet

MSstatsBioNet is an R/Bioconductor package for network analysis and
enrichment of MSstats differential abundance results in the context of
prior-knowledge biomolecular networks. It takes the output of MSstats
(or MSstatsTMT / MSstatsPTM) differential abundance analysis, maps the
analytes (proteins, PTM sites, metabolites, lipids, or drugs) to
database identifiers, queries network databases for the interactions
among them, and filters, contextualizes, and visualizes the resulting
subnetworks. Notably, it integrates with
[INDRA](https://github.com/sorgerlab/indra), a database of biological
networks assembled from the literature using text mining, enabling
interpretation of proteomic, phosphoproteomic, and metabolomic results
against past published knowledge.

MSstatsBioNet is part of the
[MSstats](https://github.com/Vitek-Lab/MSstats) family of packages,
developed and maintained by the [Vitek
Lab](https://olga-vitek-lab.khoury.northeastern.edu/) at Northeastern
University. The package and its documentation are also available at
[msstats.org](http://msstats.org).

## Installation

``` r

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("MSstatsBioNet")
```

The development version can be installed directly from this repository:

``` r

remotes::install_github("Vitek-Lab/MSstatsBioNet")
```

## Quick Start

``` r

library(MSstatsBioNet)

# Example MSstats differential abundance results (groupComparison output)
input <- data.table::fread(system.file("extdata/groupComparisonModel.csv",
                                        package = "MSstatsBioNet"))

# Describe the analytes: here UniProt accessions of proteins. For other
# molecules, change entity_type and id_type, e.g. entity_type = "metabolite"
# and id_type = "chemical_name" for compound names
entities <- prepare_entities(input, entity_type = "protein",
                             id_type = "uniprot")

# Map the identifiers to INDRA's namespaces and select significant analytes
indra <- indra_backend()
entities <- convert_ids(indra, entities)
entities <- select_entities(entities, pvalue_cutoff = 0.05)

# Retrieve the subnetwork of interactions among the selected analytes
network <- get_network(indra, entities, subnetwork_query())

head(network$nodes)
head(network$edges)

# Visualize the network (e.g. in Cytoscape, or export to HTML)
cytoscapeNetwork(network$nodes, network$edges)
```

## Input Formats

MSstatsBioNet works on **differential abundance results**, not raw
search-tool output. It accepts the `ComparisonResult` table produced by
the group-comparison functions across the MSstats ecosystem, or any
table with one row per analyte in the same format:

| Upstream package | Function producing input |
|----|----|
| [MSstats](https://github.com/Vitek-Lab/MSstats) | [`groupComparison()`](https://rdrr.io/pkg/MSstats/man/groupComparison.html) |
| [MSstatsTMT](https://github.com/Vitek-Lab/MSstatsTMT) | `groupComparisonTMT()` |
| [MSstatsPTM](https://github.com/Vitek-Lab/MSstatsPTM) | `groupComparisonPTM()` |

The input table provides, per analyte and comparison, the log2 fold
change, p-value, and adjusted p-value used for filtering and network
coloring. The analyte identifiers are read from the `Protein` column by
default; set `id_column` in
[`prepare_entities()`](https://vitek-lab.github.io/MSstatsBioNet/reference/prepare_entities.md)
to use another column.

[`prepare_entities()`](https://vitek-lab.github.io/MSstatsBioNet/reference/prepare_entities.md)
takes the type of each analyte (`entity_type`) and its identifier system
(`id_type`), either as one value for the whole table or as a column with
one value per row, so a single table can mix molecule types.
[`convert_ids()`](https://vitek-lab.github.io/MSstatsBioNet/reference/convert_ids.md)
then maps the identifiers to the backend’s namespaces. With the INDRA
backend, the supported combinations are:

| `entity_type` | `id_type` | Mapped to |
|----|----|----|
| `"protein"`, `"ptm_site"` | `"uniprot"`, `"uniprot_mnemonic"`, `"hgnc_symbol"` | HGNC |
| `"metabolite"`, `"lipid"`, `"drug"` | `"chemical_name"` | CHEBI, MESH, PUBCHEM, … (grounded by name with [Gilda](https://github.com/gyorilab/gilda)) |

`backend_capabilities(indra_backend())` lists these from R. Entities
outside the input, such as unmeasured enzymes or receptors, can be added
to a query with `get_network(include_entities = )`.

- **Databases supported:** INDRA
- **Filtering options:** p-value, fold-change, and direction filters
  ([`select_entities()`](https://vitek-lab.github.io/MSstatsBioNet/reference/select_entities.md)),
  context/topic-based filtering
  ([`filterSubnetworkByContext()`](https://vitek-lab.github.io/MSstatsBioNet/reference/filterSubnetworkByContext.md))
- **Visualization options:** Cytoscape Desktop
  ([`cytoscapeNetwork()`](https://vitek-lab.github.io/MSstatsBioNet/reference/cytoscapeNetwork.md)),
  in-browser preview
  ([`previewNetworkInBrowser()`](https://vitek-lab.github.io/MSstatsBioNet/reference/previewNetworkInBrowser.md)),
  standalone HTML export
  ([`exportNetworkToHTML()`](https://vitek-lab.github.io/MSstatsBioNet/reference/exportNetworkToHTML.md)),
  and Shiny integration
  ([`cytoscapeNetworkOutput()`](https://vitek-lab.github.io/MSstatsBioNet/reference/cytoscapeNetworkOutput.md)
  /
  [`renderCytoscapeNetwork()`](https://vitek-lab.github.io/MSstatsBioNet/reference/renderCytoscapeNetwork.md))

## Documentation

- [MSstatsBioNet
  overview](https://vitek-lab.github.io/MSstatsBioNet/vignettes/MSstatsBioNet.Rmd)
  — getting started
- [Cytoscape
  visualization](https://vitek-lab.github.io/MSstatsBioNet/vignettes/Cytoscape-Visualization.Rmd)
- [Filter by
  context](https://vitek-lab.github.io/MSstatsBioNet/vignettes/Filter-By-Context.Rmd)
- [PTM
  analysis](https://vitek-lab.github.io/MSstatsBioNet/vignettes/PTM-Analysis.Rmd)
- [Metabolomics
  analysis](https://vitek-lab.github.io/MSstatsBioNet/vignettes/Metabolomics-Analysis.Rmd)
- [Official website: msstats.org](http://msstats.org)
- [Bioconductor package page and reference
  manual](https://bioconductor.org/packages/MSstatsBioNet)

## Getting Help / Reporting Bugs

- **Questions about usage, statistical methods, or troubleshooting:**
  please post to the [MSstats Google
  Group](https://groups.google.com/forum/#!forum/msstats). This is
  monitored by the development team and searchable, so it’s the fastest
  way to get help and to see if your question has already been answered.
- **Bug reports and feature requests for this repository:** please open
  a [GitHub issue](https://github.com/Vitek-Lab/MSstatsBioNet/issues).

## References

If you use MSstatsBioNet, please cite:

1.  Wu A, Kohler D, Navada P, Robbins J, Boyle G, Boshart A, Karis K,
    Neefjes J, Konvalinka A, Sarthy J, Pino L, Gyori B, Vitek O.
    **MSstatsBioNet: Integrating Statistical Analyses with Prior
    Knowledge Biomolecular Networks for Quantitative Proteomics and
    Phosphoproteomics.** *bioRxiv*. 2026. [DOI:
    10.64898/2026.07.09.737605](https://doi.org/10.64898/2026.07.09.737605)

## Funding

MSstats development has been supported by the Chan Zuckerberg
Initiative’s [Essential Open Source Software for
Science](https://chanzuckerberg.com/eoss/proposals/).

## License

MSstatsBioNet is released under the
[Artistic-2.0](https://opensource.org/licenses/Artistic-2.0) license.
However, its dependencies may have different licenses. Notably, INDRA is
distributed under the [BSD
2-Clause](https://opensource.org/license/bsd-2-clause) license, and
INDRA’s knowledge sources may have different licenses for commercial
applications. Please refer to the [INDRA
README](https://github.com/sorgerlab/indra?tab=readme-ov-file#indra-modules)
for more information on its knowledge sources and their associated
licenses.
