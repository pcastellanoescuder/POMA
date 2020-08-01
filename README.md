
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
[![License: GPL
v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

<!-- badges: end -->

[**POMA**](http://pcastellanoescuder.github.io/POMA/) introduces a
structured, reproducible and easy use workflow for the visualization,
preprocessement, exploratory and statistical analysis of mass
spectrometry data. The main aim of `POMA` is to enable a flexible data
cleaning and statistical analysis processes in one comprehensible and
easy use R package. This package re-uses the standardized
[**MSnbase**](http://lgatto.github.io/MSnbase/) data structures,
developed by [Laurent Gatto](http://lgatto.github.io/), to achieve the
maximum flexibility and reproducibility and makes `POMA` compatible with
pre-existing [Bioconductor](https://bioconductor.org) packages.

`POMA` also offers two different Shiny app modules both for Exploratory
Data Analysis and Statistical Analysis that implements all `POMA`
functions in a easy use web interface.

  - **POMA Shiny**: Shiny version of this package.
    <https://github.com/pcastellanoescuder/POMAShiny>  
  - **POMAcounts**: Shiny version for mass spectrometry spectral counts
    data based on [Bioconductor](https://bioconductor.org) packages
    [msmsEDA](https://bioconductor.org/packages/release/bioc/html/msmsEDA.html)
    and
    [msmsTests](https://bioconductor.org/packages/release/bioc/html/msmsTests.html).
    <https://github.com/pcastellanoescuder/POMAcounts>

The [github page](https://github.com/pcastellanoescuder/POMA) is for
active development, issue tracking and forking/pulling purposes. To get
an overview of the package, see the
[*POMA-demo*](https://pcastellanoescuder.github.io/POMA/articles/POMA-demo.html)
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
BiocManager::install("pcastellanoescuder/POMA")
```

or

``` r
# install.packages("devtools")
devtools::install_github("pcastellanoescuder/POMA")
```

## Code of Conduct

Please note that the ‘POMA’ project is released with a [Contributor Code
of
Conduct](https://pcastellanoescuder.github.io/POMA/CODE_OF_CONDUCT.html).
By contributing to this project, you agree to abide by its terms.
