context("PomaLasso")

test_that("PomaLasso works", {

  data("st000336")

  imputed <- POMA::PomaImpute(st000336, method = "knn")
  normalized <- POMA::PomaNorm(imputed, method = "log_scaling")

  normalized_test <- POMA::PomaNorm(imputed, method = "log_scaling")
  Biobase::pData(normalized_test)$Group <- c(rep("C", 30), rep("G", 27))

  lasso_res <- PomaLasso(normalized, method = "lasso")
  ridge_res <- PomaLasso(normalized, method = "ridge")

  ##

  expect_error(PomaLasso(normalized, method = "lass"))
  expect_error(PomaLasso(normalized_test, method = "lasso"))
  expect_warning(PomaLasso(normalized))

  ##

  expect_false(nrow(lasso_res$coefficients) == nrow(ridge_res$coefficients))
  expect_equal(ncol(lasso_res$coefficients), ncol(ridge_res$coefficients))

  ##

  df_a <- layer_data(lasso_res$coefficientPlot)
  df_b <- layer_data(ridge_res$coefficientPlot)

  df_c <- layer_data(lasso_res$cvLassoPlot)
  df_d <- layer_data(ridge_res$cvLassoPlot)

  expect_false(length(df_a$y) == length(df_b$y))
  expect_false(length(df_c$y) == length(df_d$y))

  ##
  
  expect_error(PomaLasso(method = "lasso"))
  expect_error(PomaLasso(iris, method = "lasso"))
  
})

