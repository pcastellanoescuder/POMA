
test_that("PomaImpute handles valid SummarizedExperiment objects", {
  data <- create_mock_summarized_experiment()
  imputed_data <- PomaImpute(data, method = "mean")
  expect_is(imputed_data, "SummarizedExperiment")
})

test_that("PomaImpute stops with non-SummarizedExperiment objects", {
  data <- data.frame(matrix(runif(100), ncol = 10))
  expect_error(PomaImpute(data), "data is not a SummarizedExperiment object")
})

test_that("PomaImpute handles zeros_as_na parameter correctly", {
  data <- create_mock_summarized_experiment()
  imputed_data <- PomaImpute(data, zeros_as_na = TRUE)
  expect_is(imputed_data, "SummarizedExperiment")
})

test_that("PomaImpute handles remove_na and cutoff parameters correctly", {
  data <- create_mock_summarized_experiment()
  imputed_data <- PomaImpute(data, remove_na = TRUE, cutoff = 50)
  expect_is(imputed_data, "SummarizedExperiment")
})

test_that("PomaImpute handles different imputation methods correctly", {
  data <- create_mock_summarized_experiment()
  for (method in c("none", "half_min", "median", "mean", "min", "knn")) { # "random_forest"
    imputed_data <- PomaImpute(data, method = method)
    expect_is(imputed_data, "SummarizedExperiment")
  }
})

test_that("PomaImpute stops with incorrect method argument", {
  data <- create_mock_summarized_experiment()
  expect_error(PomaImpute(data, method = "invalid_method"), "Incorrect value for method argument")
})

