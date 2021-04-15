context("PomaLimma")

test_that("PomaLimma works", {

  data("st000284")
  data("st000336")
  
  a <- PomaLimma(st000284, contrast = "CRC-Polyp", covariates = FALSE, adjust = "fdr")
  b <- PomaLimma(st000284, contrast = "CRC-Polyp", covariates = FALSE, adjust = "bonferroni")
  c <- PomaLimma(st000284, contrast = "CRC-Polyp", covariates = TRUE, adjust = "fdr")
  d <- PomaLimma(st000284, contrast = "CRC-Polyp", covariates = TRUE, adjust = "bonferroni")

  e <- PomaLimma(st000336, contrast = "Controls-DMD", covariates = FALSE, adjust = "fdr")
  f <- PomaLimma(st000336, contrast = "Controls-DMD", covariates = TRUE, adjust = "fdr")
  
  ####

  expect_error(PomaLimma(st000284, covariates = FALSE, adjust = "fdr"))
  expect_error(PomaLimma(st000284, contrast = NULL))
  expect_warning(PomaLimma(st000284, contrast = "CRC-Polyp", covariates = FALSE))
  expect_warning(PomaLimma(st000284, contrast = "CRC-Polyp", covariates = TRUE))


  ####

  expect_equal(dim(a), dim(b))
  expect_equal(dim(b), dim(c))
  expect_equal(dim(c), dim(d))

  expect_false(all(a == b))
  expect_false(all(a == c))

  expect_false(all(b == d))

  expect_equal(dim(e), dim(f))
  
  ####

  MSnbase::pData(st000284) <- MSnbase::pData(st000284)[1]
  expect_error(PomaLimma(st000284, contrast = "CRC-Polyp", covariates = TRUE, adjust = "fdr"))
  expect_error(PomaLimma(st000284, contrast = "CRC-Polyp", covariates = TRUE, adjust = "fd"))
  
  ##
  
  expect_error(PomaLimma(contrast = "CRC-Polyp"))
  expect_error(PomaLimma(iris, contrast = "CRC-Polyp"))
  
})

