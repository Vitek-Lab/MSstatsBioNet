# MSstatsBioNet

[![Codecov test
coverage](https://codecov.io/github/Vitek-Lab/MSstatsBioNet/graph/badge.svg?token=SCPSPMTOEF)](https://codecov.io/github/Vitek-Lab/MSstatsBioNet)

This package provides a suite of functions to query various network
databases, filter queries & results, and visualize networks.

## Installation Instructions

To install this package on bioconductor, run the following command:

    if (!require("BiocManager", quietly = TRUE))
        install.packages("BiocManager")

    BiocManager::install("MSstatsBioNet")

You can install the development version of this package through Github:

    devtools::install_github("Vitek-Lab/MSstatsBioNet", build_vignettes = TRUE)

## License

This package is distributed under the
[Artistic-2.0](https://opensource.org/licenses/Artistic-2.0) license.
However, its dependencies may have different licenses.

Notably, INDRA is distributed under the [BSD
2-Clause](https://opensource.org/license/bsd-2-clause) license.
Furthermore, INDRA’s knowledge sources may have different licenses for
commercial applications. Please refer to the [INDRA
README](https://github.com/sorgerlab/indra?tab=readme-ov-file#indra-modules)
for more information on its knowledge sources and their associated
licenses.

## Databases Supported

- INDRA

## Filtering Options Supported

- P-Value Filter

## Visualization Options Supported

- Cytoscape Desktop
