
test_that("PomaCorr handles valid SummarizedExperiment objects", {
  data <- create_mock_summarized_experiment()
  result <- PomaCorr(data)
  expect_is(result, "list")
  expect_true(all(c("correlations", "corrplot") %in% names(result)))
})

test_that("PomaCorr stops with non-SummarizedExperiment objects", {
  data <- data.frame(matrix(runif(100), ncol = 10))
  expect_error(PomaCorr(data), "data is not a SummarizedExperiment object")
})

test_that("PomaCorr handles different methods correctly", {
  data <- create_mock_summarized_experiment()
  for (method in c("pearson", "kendall", "spearman")) {
    result <- PomaCorr(data, method = method)
    expect_is(result, "list")
    expect_true(all(c("correlations", "corrplot") %in% names(result)))
  }
})

test_that("PomaCorr stops with incorrect method argument", {
  data <- create_mock_summarized_experiment()
  expect_error(PomaCorr(data, method = "invalid_method"), "Incorrect value for method argument")
})

