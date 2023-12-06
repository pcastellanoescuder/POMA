
test_that("PomaPCA handles valid SummarizedExperiment objects", {
  data <- create_mock_summarized_experiment()
  pca_results <- PomaPCA(data)
  expect_is(pca_results, "list")
  expect_true(all(c("factors", "factors_plot", "eigenvalues", "eigenvalues_plot", "loadings", "loadings_plot", "biplot") %in% names(pca_results)))
})

test_that("PomaPCA stops with non-SummarizedExperiment objects", {
  data <- data.frame(matrix(runif(100), ncol = 10))
  expect_error(PomaPCA(data), "data is not a SummarizedExperiment object")
})

test_that("PomaPCA handles different parameter settings correctly", {
  data <- create_mock_summarized_experiment()
  pca_results_center_scale <- PomaPCA(data, center = TRUE, scale = TRUE)
  pca_results_ncomp <- PomaPCA(data, ncomp = 2)
  pca_results_labels <- PomaPCA(data, labels = TRUE)
  pca_results_ellipse <- PomaPCA(data, ellipse = TRUE)
  expect_is(pca_results_center_scale, "list")
  expect_is(pca_results_ncomp, "list")
  expect_is(pca_results_labels, "list")
  expect_is(pca_results_ellipse, "list")
})

