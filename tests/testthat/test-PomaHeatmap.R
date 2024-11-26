
test_that("PomaHeatmap handles valid SummarizedExperiment objects", {
  data <- create_mock_summarized_experiment()
  heatmap_plot <- PomaHeatmap(data)
  expect_is(heatmap_plot, "ggplot")
})

test_that("PomaHeatmap stops with non-SummarizedExperiment objects", {
  data <- data.frame(matrix(runif(100), ncol = 10))
  expect_error(PomaHeatmap(data), "data is not a SummarizedExperiment object")
})

test_that("PomaHeatmap handles covariates correctly", {
  data <- create_mock_summarized_experiment()
  heatmap_plot_with_covs <- PomaHeatmap(data, covs = c("group"))
  expect_is(heatmap_plot_with_covs, "ggplot")
})

test_that("PomaHeatmap handles sample_names and feature_names parameters correctly", {
  data <- create_mock_summarized_experiment()
  heatmap_plot_with_names <- PomaHeatmap(data, sample_names = FALSE, feature_names = TRUE)
  expect_is(heatmap_plot_with_names, "ggplot")
})

test_that("PomaHeatmap handles show_legend parameter correctly", {
  data <- create_mock_summarized_experiment()
  heatmap_plot_no_legend <- PomaHeatmap(data, show_legend = FALSE)
  expect_is(heatmap_plot_no_legend, "ggplot")
})

