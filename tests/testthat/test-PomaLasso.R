context("PomaLasso")

test_that("PomaLasso works", {

  data("st000336")

  normalized <- st000336 %>%
    POMA::PomaImpute(method = "knn") %>%
    POMA::PomaNorm(method = "log_scaling")

  normalized_test <- normalized
  normalized_test_less <- normalized
  SummarizedExperiment::colData(normalized_test)[,1] <- c(rep("C", 30), rep("G", 20), rep("P", 7))
  SummarizedExperiment::colData(normalized_test_less)[,1] <- "Control"

  lasso_res <- PomaLasso(normalized, alpha = 1, ntest = NULL, nfolds = 3, labels = TRUE)
  ridge_res <- PomaLasso(normalized, alpha = 0, ntest = NULL, nfolds = 10)
  enet_res <- PomaLasso(normalized, alpha = 0.5, ntest = NULL, nfolds = 5, labels = TRUE)
  lasso_self_lambda <- PomaLasso(normalized, alpha = 1, ntest = NULL, nfolds = 10, lambda = seq(0.002, 20, length.out = 100))
  
  ## TABLES

  expect_false(nrow(lasso_res$coefficients) == nrow(ridge_res$coefficients))
  expect_false(nrow(lasso_res$coefficients) == nrow(enet_res$coefficients))
  
  expect_equal(ncol(lasso_res$coefficients), ncol(ridge_res$coefficients))
  expect_equal(ncol(ridge_res$coefficients), ncol(enet_res$coefficients))
  expect_equal(ncol(lasso_res$coefficients), ncol(lasso_self_lambda$coefficients))

  ## PLOTS

  df_a <- layer_data(lasso_res$coefficientPlot)
  df_b <- layer_data(ridge_res$coefficientPlot)
  df_e <- layer_data(enet_res$coefficientPlot)

  df_c <- layer_data(lasso_res$cvLassoPlot)
  df_d <- layer_data(ridge_res$cvLassoPlot)

  expect_false(length(df_a$y) == length(df_b$y))
  expect_false(length(df_c$y) == length(df_d$y))
  expect_false(length(df_a$y) == length(df_e$y))

  ## ERRORS
  
  expect_error(PomaLasso(normalized, alpha = 2))
  expect_error(PomaLasso(normalized, alpha = -0.5))
  expect_error(PomaLasso(iris, alpha = 1))
  expect_error(PomaLasso(normalized_test, alpha = 1))
  expect_error(PomaLasso(normalized_test_less, alpha = 1))
  expect_error(PomaLasso())
  expect_error(PomaLasso(normalized, ntest = 60))
  expect_error(PomaLasso(normalized, ntest = 2))
  
})

