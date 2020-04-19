context("PomaLimma")

test_that("PomaLimma works", {

  data("st000284")

  a <- PomaLimma(st000284, contrast = "C-H", covariates = FALSE, adjust = "fdr")
  b <- PomaLimma(st000284, contrast = "C-H", covariates = FALSE, adjust = "bonferroni")
  c <- PomaLimma(st000284, contrast = "C-H", covariates = TRUE, adjust = "fdr")
  d <- PomaLimma(st000284, contrast = "C-H", covariates = TRUE, adjust = "bonferroni")

  ####

  expect_error(PomaLimma(st000284, covariates = FALSE, adjust = "fdr"))
  expect_error(PomaLimma(st000284, contrast = NULL))
  expect_warning(PomaLimma(st000284, contrast = "C-H", covariates = FALSE))
  expect_warning(PomaLimma(st000284, contrast = "C-H", covariates = TRUE))


  ####

  expect_equal(dim(a), dim(b))
  expect_equal(dim(b), dim(c))
  expect_equal(dim(c), dim(d))

  expect_false(all(a == b))
  expect_false(all(a == c))

  expect_false(all(b == d))

  ####

  Biobase::pData(st000284) <- Biobase::pData(st000284)[1]
  expect_error(PomaLimma(st000284, contrast = "C-H", covariates = TRUE, adjust = "fdr"))
  expect_error(PomaLimma(st000284, contrast = "C-H", covariates = TRUE, adjust = "fd"))
  
})

