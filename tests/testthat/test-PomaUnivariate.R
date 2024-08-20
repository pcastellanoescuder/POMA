
test_that("PomaUnivariate handles valid SummarizedExperiment objects", {
  data <- create_mock_summarized_experiment(binary = TRUE)
  univariate_results <- PomaUnivariate(data, method = "ttest")
  expect_is(univariate_results, "tbl_df")
})

test_that("PomaUnivariate stops with non-SummarizedExperiment objects", {
  data <- data.frame(matrix(runif(100), ncol = 10))
  expect_error(PomaUnivariate(data), "data is not a SummarizedExperiment object")
})

test_that("PomaUnivariate handles different methods correctly for 2 groups", {
  data <- create_mock_summarized_experiment(binary = TRUE)
  for (method in c("ttest", "mann")) {
    univariate_results <- PomaUnivariate(data, method = method)
    expect_is(univariate_results, "tbl_df")
  }
})

test_that("PomaUnivariate handles different methods correctly for more than 2 groups", {
  data <- create_mock_summarized_experiment()
  for (method in c("anova", "kruskal")) {
    univariate_results <- PomaUnivariate(data, method = method)
    expect_is(univariate_results, "list")
    expect_true("result" %in% names(univariate_results))
  }
})

test_that("PomaUnivariate stops with incorrect method argument", {
  data <- create_mock_summarized_experiment()
  expect_error(PomaUnivariate(data, method = "invalid_method"), "Incorrect value for method argument")
})

test_that("PomaUnivariate handles paired, var_equal, and adjust parameters correctly", {
  data_binary <- create_mock_summarized_experiment(binary = TRUE)
  data_binary_paired <- create_mock_summarized_experiment(binary = TRUE, paired = TRUE)

  univariate_results_paired <- PomaUnivariate(data_binary_paired, method = "ttest", paired = TRUE)
  univariate_results_var_equal <- PomaUnivariate(data_binary, method = "ttest", var_equal = TRUE)
  univariate_results_adjusted <- PomaUnivariate(data_binary, method = "ttest", adjust = "holm")
  expect_is(univariate_results_paired, "tbl_df")
  expect_is(univariate_results_var_equal, "tbl_df")
  expect_is(univariate_results_adjusted, "tbl_df")
})

