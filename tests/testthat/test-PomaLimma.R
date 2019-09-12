context("PomaLimma")

test_that("PomaLimma works", {

  library(tidyverse)
  library(limma)

  data <- vroom::vroom("data_ST000284/MET_CRC_ST000284.csv", delim = ",")
  covariates <- vroom::vroom("data_ST000284/COV_CRC_ST000284.csv", delim = ",")

  a <- PomaLimma(data, contrast = "C-H", covariates = NULL, adjust = "fdr")
  b <- PomaLimma(data, contrast = "C-H", covariates = NULL, adjust = "bonferroni")
  c <- PomaLimma(data, contrast = "C-H", covariates = covariates, adjust = "fdr")
  d <- PomaLimma(data, contrast = "C-H", covariates = covariates, adjust = "bonferroni")

  ####

  expect_error(PomaLimma(data, covariates = NULL, adjust = "fdr"))
  expect_error(PomaLimma(data, contrast = NULL))
  expect_warning(PomaLimma(data, contrast = "C-H", covariates = NULL))
  expect_warning(PomaLimma(data, contrast = "C-H", covariates = covariates))


  ##########

  expect_equal(dim(a), dim(b))
  expect_equal(dim(b), dim(c))
  expect_equal(dim(c), dim(d))

  expect_false(all(a == b))
  expect_false(all(a == c))

  expect_false(all(b == d))

})

