
test_that("PomaCreateObject handles features and metadata correctly", {
  mock_data <- create_mock_data()
  se_object <- PomaCreateObject(metadata = mock_data$metadata, features = mock_data$features)
  expect_is(se_object, "SummarizedExperiment")
})

test_that("PomaCreateObject handles missing metadata", {
  mock_data <- create_mock_data()
  se_object <- PomaCreateObject(features = mock_data$features)
  expect_is(se_object, "SummarizedExperiment")
})

test_that("PomaCreateObject handles factor_levels parameter correctly", {
  mock_data <- create_mock_data()
  mock_data$metadata$numeric_var <- 1:10
  se_object <- PomaCreateObject(metadata = mock_data$metadata, features = mock_data$features, factor_levels = 5)
  expect_is(se_object, "SummarizedExperiment")
  expect_true("numeric_var" %in% colnames(SummarizedExperiment::colData(se_object)))
  expect_true(is.numeric(SummarizedExperiment::colData(se_object)$numeric_var))
})

