
test_that("PomaPLS handles valid SummarizedExperiment objects", {
  data <- create_mock_summarized_experiment()
  data$numeric_var <- 1:20
  pls_results <- PomaPLS(data, method = "pls")
  expect_is(pls_results, "list")
  expect_true(all(c("factors", "factors_plot", "loadings", "loadings_plot") %in% names(pls_results)))
})

test_that("PomaPLS stops with non-SummarizedExperiment objects", {
  data <- data.frame(matrix(runif(100), ncol = 10))
  expect_error(PomaPLS(data, method = "pls"), "data is not a SummarizedExperiment object")
})

test_that("PomaPLS handles different PLS methods correctly", {
  data <- create_mock_summarized_experiment()
  plsda_results <- PomaPLS(data, method = "plsda", ncomp = 3)
  splsda_results <- PomaPLS(data, method = "splsda", ncomp = 3, num_features = 3)
  expect_is(plsda_results, "list")
  expect_is(splsda_results, "list")
})

test_that("PomaPLS handles different parameter settings correctly", {
  data <- create_mock_summarized_experiment()
  data$numeric_var <- 1:20
  pls_results_ncomp <- PomaPLS(data, method = "pls", ncomp = 3)
  plsda_results_labels <- PomaPLS(data, method = "plsda", ncomp = 3, labels = TRUE)
  splsda_results_ellipse <- PomaPLS(data, method = "splsda", ncomp = 3, num_features = 3, ellipse = TRUE)
  expect_is(pls_results_ncomp, "list")
  expect_is(plsda_results_labels, "list")
  expect_is(splsda_results_ellipse, "list")
})

test_that("PomaPLS stops with incorrect method argument", {
  data <- create_mock_summarized_experiment()
  expect_error(PomaPLS(data, method = "invalid_method"), "Incorrect value for method argument")
})

