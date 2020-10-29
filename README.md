
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
[![Bioc
Status](https://bioconductor.org/shields/build/devel/bioc/POMA.svg)](https://bioconductor.org/checkResults/devel/bioc-LATEST/POMA/)
[![Codecov test
coverage](https://codecov.io/gh/pcastellanoescuder/POMA/branch/master/graph/badge.svg)](https://codecov.io/gh/pcastellanoescuder/POMA?branch=master)
[![Last
Commit](https://img.shields.io/github/last-commit/pcastellanoescuder/POMA.svg)](https://github.com/pcastellanoescuder/POMA/commits/master)
[![License: GPL
v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

<!-- badges: end -->

[**POMA**](http://pcastellanoescuder.github.io/POMA/) introduces a
structured, reproducible and easy-to-use workflow for the visualization,
pre-processing, exploratory and statistical analysis of mass
spectrometry data. The main aim of `POMA` is to enable a flexible data
cleaning and statistical analysis processes in one comprehensible and
user-friendly R package. This package uses the standardized
[**MSnbase**](http://lgatto.github.io/MSnbase/) data structures,
developed by [Laurent Gatto](http://lgatto.github.io/), to achieve the
maximum flexibility and reproducibility and makes `POMA` compatible with
other [Bioconductor](https://bioconductor.org) packages.

`POMA` also has two different Shiny app modules both for Exploratory
Data Analysis and Statistical Analysis that implement all `POMA`
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
