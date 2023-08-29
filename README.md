
# POMA <img src='man/figures/logo.png' align="right" height="139" />

<!-- badges: start -->

[![Lifecycle:
stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://www.tidyverse.org/lifecycle/#stable)
[![Actions
Status](https://github.com/pcastellanoescuder/POMA/workflows/R-CMD-check/badge.svg)](https://github.com/pcastellanoescuder/POMA/actions)
[![Codecov test
coverage](https://codecov.io/gh/pcastellanoescuder/POMA/branch/master/graph/badge.svg)](https://codecov.io/gh/pcastellanoescuder/POMA?branch=master)
[![CodeFactor](https://www.codefactor.io/repository/github/pcastellanoescuder/POMA/badge)](https://www.codefactor.io/repository/github/pcastellanoescuder/POMA)
[![Last
Commit](https://img.shields.io/github/last-commit/pcastellanoescuder/POMA.svg)](https://github.com/pcastellanoescuder/POMA/commits/master)
[![License: GPL
v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

| *BioC* branch                                                           | Status                                                                                                                                                  | Version                                                                                                                                            | Dependencies                                                                                                                                         | Rank                                                                                                                         |
|-------------------------------------------------------------------------|---------------------------------------------------------------------------------------------------------------------------------------------------------|----------------------------------------------------------------------------------------------------------------------------------------------------|------------------------------------------------------------------------------------------------------------------------------------------------------|------------------------------------------------------------------------------------------------------------------------------|
| [Release](http://bioconductor.org/packages/release/bioc/html/POMA.html) | [![Bioc release status](https://bioconductor.org/shields/build/release/bioc/POMA.svg)](https://bioconductor.org/checkResults/release/bioc-LATEST/POMA/) | [![BioC released version](https://img.shields.io/badge/release%20version-1.6.0-blue.svg)](https://www.bioconductor.org/packages/POMA)              | [![Dependencies](http://bioconductor.org/shields/dependencies/release/POMA.svg)](http://bioconductor.org/packages/release/bioc/html/POMA.html#since) | [![Rank](http://www.bioconductor.org/shields/downloads/release/POMA.svg)](https://bioconductor.org/packages/stats/bioc/POMA) |
| [Devel](http://bioconductor.org/packages/devel/bioc/html/POMA.html)     | [![Bioc devel status](https://bioconductor.org/shields/build/devel/bioc/POMA.svg)](https://bioconductor.org/checkResults/devel/bioc-LATEST/POMA/)       | [![BioC devel version](https://img.shields.io/badge/devel%20version-1.8.27-blue.svg)](https://bioconductor.org/packages/devel/bioc/html/POMA.html) | [![Dependencies](http://bioconductor.org/shields/dependencies/devel/POMA.svg)](http://bioconductor.org/packages/devel/bioc/html/POMA.html#since)     | [![Rank](http://www.bioconductor.org/shields/downloads/devel/POMA.svg)](https://bioconductor.org/packages/stats/bioc/POMA)   |

<!-- badges: end -->

[**POMA**](http://pcastellanoescuder.github.io/POMA/) introduces a
reproducible and easy-to-use toolkit for visualization, pre-processing,
exploration, and statistical analysis of omics datasets. The main aim of
POMA is to enable a flexible data cleaning and statistical analysis
processes in one comprehensible and user-friendly R package. This
package uses the standardized
[**SummarizedExperiment**](https://bioconductor.org/packages/release/bioc/html/SummarizedExperiment.html)
class to achieve the maximum flexibility and reproducibility with other
[Bioconductor](https://bioconductor.org) packages.

POMA provides two different Shiny apps both for exploratory data
analysis and statistical analysis that implement all POMA functions in
two user-friendly web interfaces.

- **POMAShiny**: Shiny version of this package.
  <https://github.com/pcastellanoescuder/POMAShiny>  
- **POMAcounts**: Shiny version for exploratory and statistical analysis
  of mass spectrometry spectral counts data and RNAseq data.
  <https://github.com/pcastellanoescuder/POMAcounts>

The [GitHub page](https://github.com/pcastellanoescuder/POMA) is for
active development, issue tracking and forking/pulling purposes. To get
an overview of the package, see the [*POMA
Workflow*](https://pcastellanoescuder.github.io/POMA/articles/POMA-demo.html)
vignette.

## Installation

To install Bioconductor version:

``` r
# install.packages("BiocManager")
BiocManager::install("POMA")
```

If you need the GitHub version (not recommended), use:

``` r
# install.packages("devtools")
devtools::install_github("pcastellanoescuder/POMA")
```

## Citation

Castellano-Escuder P, González-Domínguez R, Carmona-Pontaque F, et
al. POMAShiny: A user-friendly web-based workflow for metabolomics and
proteomics data analysis. PLoS Comput Biol. 2021 Jul 1;17(7):e1009148.
doi: 10.1371/journal.pcbi.1009148. PMID: 34197462; PMCID: PMC8279420.

### Cited In

Bellio C, Emperador M, Castellano P, et al. GDF15 Is an Eribulin
Response Biomarker also Required for Survival of DTP Breast Cancer
Cells. Cancers (Basel). 2022 May 23;14(10):2562. doi:
10.3390/cancers14102562. PMID: 35626166; PMCID: PMC9139899.

González-Domínguez R, Castellano-Escuder P, Lefèvre-Arbogast S, et
al. Apolipoprotein E and sex modulate fatty acid metabolism in a
prospective observational study of cognitive decline. Alzheimers Res
Ther. 2022 Jan 3;14(1):1. doi: 10.1186/s13195-021-00948-8. PMID:
34980257; PMCID: PMC8725342.

González-Domínguez R, Castellano-Escuder P, Carmona F, et al. Food and
Microbiota Metabolites Associate with Cognitive Decline in Older
Subjects: A 12-Year Prospective Study. Mol Nutr Food Res. 2021
Dec;65(23):e2100606. doi: 10.1002/mnfr.202100606. Epub 2021 Oct 28.
PMID: 34661340.

Peron G, Gargari G, Meroño T, et al. Crosstalk among intestinal barrier,
gut microbiota and serum metabolome after a polyphenol-rich diet in
older subjects with “leaky gut”: The MaPLE trial. Clin Nutr. 2021
Oct;40(10):5288-5297. doi: 10.1016/j.clnu.2021.08.027. Epub 2021 Sep 9.
PMID: 34534897.

## News

Click
[here](https://github.com/pcastellanoescuder/POMA/blob/master/NEWS.md)
for the latest package updates.

## Code of Conduct

Please note that the ‘POMA’ project is released with a [Contributor Code
of
Conduct](https://pcastellanoescuder.github.io/POMA/CODE_OF_CONDUCT.html).
By contributing to this project, you agree to abide by its terms.
