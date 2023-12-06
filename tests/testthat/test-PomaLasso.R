
test_that("PomaLasso handles valid SummarizedExperiment objects", {
  data <- create_mock_summarized_experiment(binary = TRUE)
  lasso_results <- PomaLasso(data, alpha = 1)
  expect_is(lasso_results, "list")
  expect_true(all(c("coefficients", "coefficients_plot", "cv_plot", "model") %in% names(lasso_results)))
})

test_that("PomaLasso stops with non-SummarizedExperiment objects", {
  data <- data.frame(matrix(runif(100), ncol = 10))
  expect_error(PomaLasso(data), "data is not a SummarizedExperiment object")
})

test_that("PomaLasso handles different alpha values correctly", {
  data <- create_mock_summarized_experiment(binary = TRUE)
  for (alpha in c(0, 0.5, 1)) { # Testing Ridge, Elasticnet, and Lasso
    lasso_results <- PomaLasso(data, alpha = alpha)
    expect_is(lasso_results, "list")
    expect_true(all(c("coefficients", "coefficients_plot", "cv_plot", "model") %in% names(lasso_results)))
  }
})

test_that("PomaLasso stops with incorrect alpha argument", {
  data <- create_mock_summarized_experiment(binary = TRUE)
  expect_error(PomaLasso(data, alpha = -1), "alpha must be a number between 0 and 1")
  expect_error(PomaLasso(data, alpha = 2), "alpha must be a number between 0 and 1")
})

test_that("PomaLasso handles ntest and nfolds parameters correctly", {
  data <- create_mock_summarized_experiment(binary = TRUE)
  lasso_results_with_ntest <- PomaLasso(data, ntest = 10)
  lasso_results_with_nfolds <- PomaLasso(data, nfolds = 5)
  expect_is(lasso_results_with_ntest, "list")
  expect_true("confusion_matrix" %in% names(lasso_results_with_ntest))
  expect_is(lasso_results_with_nfolds, "list")
})

