
test_that("PomaNorm handles valid SummarizedExperiment objects", {
  data <- create_mock_summarized_experiment()
  normalized_data <- PomaNorm(data, method = "auto_scaling")
  expect_is(normalized_data, "SummarizedExperiment")
})

test_that("PomaNorm stops with non-SummarizedExperiment objects", {
  data <- data.frame(matrix(runif(100), ncol = 10))
  expect_error(PomaNorm(data), "data is not a SummarizedExperiment object")
})

test_that("PomaNorm handles sample normalization methods correctly", {
  data <- create_mock_summarized_experiment()
  normalized_data_sum <- PomaNorm(data, sample_norm = "sum")
  normalized_data_quantile <- PomaNorm(data, sample_norm = "quantile")
  expect_is(normalized_data_sum, "SummarizedExperiment")
  expect_is(normalized_data_quantile, "SummarizedExperiment")
})

test_that("PomaNorm handles different normalization methods correctly", {
  data <- create_mock_summarized_experiment()
  for (method in c("none", "auto_scaling", "level_scaling", "log_scaling", "log",
                   "vast_scaling", "log_pareto", "min_max", "box_cox")) {
    normalized_data <- PomaNorm(data, method = method)
    expect_is(normalized_data, "SummarizedExperiment")
  }
})

test_that("PomaNorm stops with incorrect normalization method argument", {
  data <- create_mock_summarized_experiment()
  expect_error(PomaNorm(data, method = "invalid_method"), "Incorrect value for method argument")
})

