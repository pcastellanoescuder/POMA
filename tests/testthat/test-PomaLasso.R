context("PomaLasso")

test_that("PomaLasso works", {

  data("st000284")

  imputed <- POMA::PomaImpute(st000284, method = "knn")
  normalized <- POMA::PomaNorm(imputed, method = "log_scaling")
  # normalized_test <- POMA::PomaNorm(imputed, method = "log_scaling")
  # Biobase::pData(normalized_test)$group <- c(rep("C", 100), rep("H", 29), rep("G", 3))

  ##

  expect_error(PomaLasso(normalized, method = "lass"))
  # expect_error(PomaLasso(normalized_test, method = "lasso"))
  expect_warning(PomaLasso(normalized))

})

