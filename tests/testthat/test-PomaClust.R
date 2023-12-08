
test_that("PomaClust handles valid SummarizedExperiment objects", {
  data <- create_mock_summarized_experiment()
  result <- PomaClust(data)
  expect_is(result, "list")
  expect_true(all(c("mds_coordinates", "mds_plot", "optimal_clusters_number", "optimal_clusters_plot") %in% names(result)))
})

test_that("PomaClust stops with non-SummarizedExperiment objects", {
  data <- data.frame(matrix(runif(100), ncol = 10))
  expect_error(PomaClust(data), "data is not a SummarizedExperiment object")
})

test_that("PomaClust handles different methods correctly", {
  data <- create_mock_summarized_experiment()
  for (method in c("euclidean", "maximum", "manhattan", "canberra", "minkowski")) {
    result <- PomaClust(data, method = method)
    expect_is(result, "list")
  }
})

test_that("PomaClust stops with incorrect method argument", {
  data <- create_mock_summarized_experiment()
  expect_error(PomaClust(data, method = "invalid_method"), "Incorrect value for method argument")
})

test_that("PomaClust handles show_clusters and labels parameters correctly", {
  data <- create_mock_summarized_experiment()
  result_with_clusters <- PomaClust(data, show_clusters = TRUE)
  result_with_labels <- PomaClust(data, labels = TRUE)
  expect_is(result_with_clusters, "list")
  expect_is(result_with_labels, "list")
})

