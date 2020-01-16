
# POMA <img src='man/figures/logo.png' align="right" height="139" />

<!-- badges: start -->

[![Lifecycle:
stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://www.tidyverse.org/lifecycle/#stable)
[![Build
Status](https://travis-ci.org/pcastellanoescuder/POMA.svg?branch=master)](https://travis-ci.org/pcastellanoescuder/POMA)
[![AppVeyor build
status](https://ci.appveyor.com/api/projects/status/github/pcastellanoescuder/POMA?branch=master&svg=true)](https://ci.appveyor.com/project/pcastellanoescuder/POMA)
[![Codecov test
coverage](https://codecov.io/gh/pcastellanoescuder/POMA/branch/master/graph/badge.svg)](https://codecov.io/gh/pcastellanoescuder/POMA?branch=master)
[![License: GPL
v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

<!-- badges: end -->

[**POMA**](http://pcastellanoescuder.github.io/POMA/) is an R package
that contains different robust functions for processing, exploratory
analysis and statistical analysis of mass spectrometry data, such as
metabolomics or proteomics. This package is based on the standardized
[**MSnbase**](http://lgatto.github.io/MSnbase/) infrastructure,
developed by [Laurent Gatto](http://lgatto.github.io/) in October 2010,
to achieve the maximum flexibility and reproducibility and makes `POMA`
compatible with pre-existing [Bioconductor](https://bioconductor.org)
packages .

`POMA` also offers two different Shiny app modules
(<http://poma.netlify.com>) both for Exploratory Data Analysis and
Statistical Analysis that implements all `POMA` functions in an easy use
web interface.

The [github page](https://github.com/pcastellanoescuder/POMA) page is
for active development, issue tracking and forking/pulling purposes. To
get an overview of the package, see the
[*POMA-demo*](https://pcastellanoescuder.github.io/POMA/articles/POMA-demo.html)
vignette.

## Installation

``` r
# install.packages("devtools")
devtools::install_github("pcastellanoescuder/POMA")
```

## Code of Conduct

Please note that the ‘POMA’ project is released with a [Contributor Code
of
Conduct](https://github.com/pcastellanoescuder/POMA/blob/master/CODE_OF_CONDUCT.md).
By contributing to this project, you agree to abide by its terms.
