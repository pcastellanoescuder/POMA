context("PomaUnivariate")

test_that("PomaUnivariate works", {

  # library(tidyverse)

  data("st000284")

  dims_for_ttest <- c(ncol(t(Biobase::exprs(st000284))), 6)
  dims_for_aov <- c(ncol(t(Biobase::exprs(st000284))) + ncol(Biobase::pData(st000284)) -1, 2)

  univ_ttest <- PomaUnivariate(st000284, method = "ttest", adjust = "fdr")
  univ_aov_cov <- PomaUnivariate(st000284, covariates = T, method = "anova", adjust = "fdr")
  univ_aov_cov2 <- PomaUnivariate(st000284, covariates = T, method = "anova", adjust = "bonferroni")

  univ_mann <- PomaUnivariate(st000284, method = "mann", adjust = "fdr")
  univ_kruskal <- PomaUnivariate(st000284, method = "kruskal", adjust = "fdr")
  univ_kruskal2 <- PomaUnivariate(st000284, method = "kruskal", adjust = "BY")

  ####

  expect_equal(dims_for_ttest, dim(univ_ttest))
  expect_equal(dims_for_aov, dim(univ_aov_cov))
  expect_false(all(univ_aov_cov == univ_aov_cov2))

  expect_error(PomaUnivariate(st000284, covariates = T, method = "anov", adjust = "fdr"))
  expect_warning(PomaUnivariate(st000284, covariates = T, method = "anova"))
  expect_error(PomaUnivariate(st000284, covariates = T, adjust = "fdr"))

  expect_equal(dims_for_ttest, dim(univ_mann))

  expect_equal(dim(univ_kruskal), dim(univ_kruskal2))
  expect_false(all(univ_kruskal == univ_kruskal2))

})

