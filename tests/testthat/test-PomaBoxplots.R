
test_that("PomaBoxplots handles valid SummarizedExperiment objects", {
  data <- create_mock_summarized_experiment()
  plot_samples <- PomaBoxplots(data, x = "samples")
  plot_features <- PomaBoxplots(data, x = "features")
  expect_is(plot_samples, "ggplot")
  expect_is(plot_features, "ggplot")
})

test_that("PomaBoxplots stops with non-SummarizedExperiment objects", {
  data <- data.frame(matrix(runif(100), ncol = 10))
  expect_error(PomaBoxplots(data), "data is not a SummarizedExperiment object")
})

test_that("PomaBoxplots handles violin plot option correctly", {
  data <- create_mock_summarized_experiment()
  plot_violin <- PomaBoxplots(data, violin = TRUE)
  expect_is(plot_violin, "ggplot")
})

test_that("PomaBoxplots stops with incorrect x argument", {
  data <- create_mock_summarized_experiment()
  expect_error(PomaBoxplots(data, x = "invalid_option"), "Incorrect value for x argument")
})

test_that("PomaBoxplots handles feature_name parameter correctly", {
  data <- create_mock_summarized_experiment()
  plot_specific_feature <- PomaBoxplots(data, x = "features", feature_name = c("V2"))
  expect_is(plot_specific_feature, "ggplot")
})

test_that("PomaBoxplots stops with non-existing feature names", {
  data <- create_mock_summarized_experiment()
  expect_error(PomaBoxplots(data, x = "features", feature_name = c("non_existing_feature")), "Features not found")
})

test_that("PomaBoxplots applies theme parameters correctly", {
  data <- create_mock_summarized_experiment()
  plot_with_theme <- PomaBoxplots(data, theme_params = list(legend_title = TRUE))
  expect_is(plot_with_theme, "ggplot")
})

