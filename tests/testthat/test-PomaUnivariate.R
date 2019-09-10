context("PomaUnivariate")

test_that("PomaUnivariate works", {

  library(tidyverse)

  data <- vroom::vroom("data_ST000284/MET_CRC_ST000284.csv", delim = ",")
  covariates <- vroom::vroom("data_ST000284/COV_CRC_ST000284.csv", delim = ",")

  dims_for_ttest <- c((ncol(data)-2), 6)
  dims_for_aov <- c((ncol(data)-2), (2*(ncol(covariates)-1)) +2)

  univ_ttest <- PomaUnivariate(data, method = "ttest", adjust = "fdr")
  univ_aov_cov <- PomaUnivariate(data, covariates, method = "anova", adjust = "fdr")
  univ_aov_cov2 <- PomaUnivariate(data, covariates, method = "anova", adjust = "bonferroni")

  univ_mann <- PomaUnivariate(data, method = "mann", adjust = "fdr")
  univ_kruskal <- PomaUnivariate(data, method = "kruskal", adjust = "fdr")
  univ_kruskal2 <- PomaUnivariate(data, method = "kruskal", adjust = "BY")

  ####

  expect_equal(dims_for_ttest, dim(univ_ttest))
  expect_equal(dims_for_aov, dim(univ_aov_cov))
  expect_false(all(univ_aov_cov == univ_aov_cov2))

  expect_error(PomaUnivariate(data, covariates, method = "anov", adjust = "fdr"))
  expect_warning(PomaUnivariate(data, covariates, method = "anova"))
  expect_error(PomaUnivariate(data, covariates, adjust = "fdr"))

  expect_equal(dims_for_ttest, dim(univ_mann))

  expect_equal(dim(univ_kruskal), dim(univ_kruskal2))
  expect_false(all(univ_kruskal == univ_kruskal2))

})

