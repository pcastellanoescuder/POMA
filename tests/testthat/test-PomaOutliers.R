
test_that("PomaOutliers handles valid SummarizedExperiment objects", {
  data <- create_mock_summarized_experiment()
  result <- PomaOutliers(data)
  expect_is(result, "list")
  expect_true(all(c("polygon_plot", "distance_boxplot", "outliers", "data") %in% names(result)))
})

test_that("PomaOutliers stops with non-SummarizedExperiment objects", {
  data <- data.frame(matrix(runif(100), ncol = 10))
  expect_error(PomaOutliers(data), "data is not a SummarizedExperiment object")
})

test_that("PomaOutliers stops with invalid method argument", {
  data <- create_mock_summarized_experiment()
  expect_error(PomaOutliers(data, method = "invalid_method"), "Incorrect value for method argument")
})

test_that("PomaOutliers stops with invalid type argument", {
  data <- create_mock_summarized_experiment()
  expect_error(PomaOutliers(data, type = "invalid_type"), "Incorrect value for type argument")
})

test_that("PomaOutliers returns expected results with different methods", {
  data <- create_mock_summarized_experiment()
  for (method in c("euclidean", "maximum", "manhattan", "canberra", "minkowski")) {
    result <- PomaOutliers(data, method = method)
    expect_is(result, "list")
  }
})

test_that("PomaOutliers handles different outlier coefficients", {
  data <- create_mock_summarized_experiment()
  for (coef in seq(1, 5, by = 1)) {
    result <- PomaOutliers(data, coef = coef)
    expect_is(result, "list")
  }
})

test_that("PomaOutliers handles labels parameter correctly", {
  data <- create_mock_summarized_experiment()
  result_with_labels <- PomaOutliers(data, labels = TRUE)
  result_without_labels <- PomaOutliers(data, labels = FALSE)
  expect_is(result_with_labels, "list")
  expect_is(result_without_labels, "list")
})

