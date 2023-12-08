
test_that("PomaDensity handles valid SummarizedExperiment objects", {
  data <- create_mock_summarized_experiment()
  plot_samples <- PomaDensity(data, x = "samples")
  plot_features <- PomaDensity(data, x = "features")
  expect_is(plot_samples, "ggplot")
  expect_is(plot_features, "ggplot")
})

test_that("PomaDensity stops with non-SummarizedExperiment objects", {
  data <- data.frame(matrix(runif(100), ncol = 10))
  expect_error(PomaDensity(data), "data is not a SummarizedExperiment object")
})

test_that("PomaDensity stops with incorrect x argument", {
  data <- create_mock_summarized_experiment()
  expect_error(PomaDensity(data, x = "invalid_option"), "Incorrect value for x argument")
})

test_that("PomaDensity handles feature_name parameter correctly", {
  data <- create_mock_summarized_experiment()
  plot_specific_feature <- PomaDensity(data, x = "features", feature_name = c("V1"))
  expect_is(plot_specific_feature, "ggplot")
})

test_that("PomaDensity stops with non-existing feature names", {
  data <- create_mock_summarized_experiment()
  expect_error(PomaDensity(data, x = "features", feature_name = c("non_existing_feature")), "Features not found")
})

test_that("PomaDensity applies theme parameters correctly", {
  data <- create_mock_summarized_experiment()
  plot_with_theme <- PomaDensity(data, theme_params = list(legend_title = TRUE))
  expect_is(plot_with_theme, "ggplot")
})

