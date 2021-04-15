
# POMA <img src='man/figures/logo.png' align="right" height="139" />

<!-- badges: start -->

[![Lifecycle:
stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://www.tidyverse.org/lifecycle/#stable)
[![Build
Status](https://travis-ci.org/pcastellanoescuder/POMA.svg?branch=master)](https://travis-ci.org/pcastellanoescuder/POMA)
[![AppVeyor build
status](https://ci.appveyor.com/api/projects/status/github/pcastellanoescuder/POMA?branch=master&svg=true)](https://ci.appveyor.com/project/pcastellanoescuder/POMA)
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
| ----------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------- |
| [Release](http://bioconductor.org/packages/release/bioc/html/POMA.html) | [![Bioc release status](https://bioconductor.org/shields/build/release/bioc/POMA.svg)](https://bioconductor.org/checkResults/release/bioc-LATEST/POMA/) | [![BioC released version](https://img.shields.io/badge/release%20version-1.0.0-blue.svg)](https://www.bioconductor.org/packages/POMA)              | [![Dependencies](http://bioconductor.org/shields/dependencies/release/POMA.svg)](http://bioconductor.org/packages/release/bioc/html/POMA.html#since) | [![Rank](http://www.bioconductor.org/shields/downloads/release/POMA.svg)](https://bioconductor.org/packages/stats/bioc/POMA) |
| [Devel](http://bioconductor.org/packages/devel/bioc/html/POMA.html)     | [![Bioc devel status](https://bioconductor.org/shields/build/devel/bioc/POMA.svg)](https://bioconductor.org/checkResults/devel/bioc-LATEST/POMA/)       | [![BioC devel version](https://img.shields.io/badge/devel%20version-1.1.11-blue.svg)](https://bioconductor.org/packages/devel/bioc/html/POMA.html) | [![Dependencies](http://bioconductor.org/shields/dependencies/devel/POMA.svg)](http://bioconductor.org/packages/devel/bioc/html/POMA.html#since)     | [![Rank](http://www.bioconductor.org/shields/downloads/devel/POMA.svg)](https://bioconductor.org/packages/stats/bioc/POMA)   |

<!-- badges: end -->

[**POMA**](http://pcastellanoescuder.github.io/POMA/) introduces a
structured, reproducible and easy-to-use workflow for the visualization,
pre-processing, exploration, and statistical analysis of mass
spectrometry data. The main aim of `POMA` is to enable a flexible data
cleaning and statistical analysis processes in one comprehensible and
user-friendly R package. This package uses the standardized
[**MSnbase**](http://lgatto.github.io/MSnbase/) data structures, to
achieve the maximum flexibility and reproducibility and makes `POMA`
compatible with other [Bioconductor](https://bioconductor.org) packages.

`POMA` also has two different Shiny app modules both for exploratory
data analysis and statistical analysis that implement all `POMA`
functions in two user-friendly web interfaces.

  - **POMAShiny**: Shiny version of this package.
    <https://github.com/pcastellanoescuder/POMAShiny>  
  - **POMAcounts**: Shiny version for mass spectrometry spectral counts
    data based on [Bioconductor](https://bioconductor.org) packages
    [msmsEDA](https://bioconductor.org/packages/release/bioc/html/msmsEDA.html)
    and
    [msmsTests](https://bioconductor.org/packages/release/bioc/html/msmsTests.html).
    <https://github.com/pcastellanoescuder/POMAcounts>

The [github page](https://github.com/pcastellanoescuder/POMA) is for
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

If you need the GitHub version (not recommended unless you know what you
are doing), use:

``` r
# install.packages("devtools")
devtools::install_github("pcastellanoescuder/POMA")
```

## Code of Conduct

Please note that the ‘POMA’ project is released with a [Contributor Code
of
Conduct](https://pcastellanoescuder.github.io/POMA/CODE_OF_CONDUCT.html).
By contributing to this project, you agree to abide by its terms.
