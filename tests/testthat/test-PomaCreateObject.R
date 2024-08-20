
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
  mock_data$metadata$numeric_var <- c(rep(1:5), rep(1:5))
  se_object <- PomaCreateObject(metadata = mock_data$metadata, features = mock_data$features, factor_levels = 6)
  expect_is(se_object, "SummarizedExperiment")
  expect_true("numeric_var" %in% colnames(SummarizedExperiment::colData(se_object)))
  expect_true(is.factor(SummarizedExperiment::colData(se_object)$numeric_var))
  
  mock_data <- create_mock_data()
  mock_data$metadata$numeric_var <- 1:10
  se_object <- PomaCreateObject(metadata = mock_data$metadata, features = mock_data$features, factor_levels = 6)
  expect_is(se_object, "SummarizedExperiment")
  expect_true("numeric_var" %in% colnames(SummarizedExperiment::colData(se_object)))
  expect_true(is.numeric(SummarizedExperiment::colData(se_object)$numeric_var))
})

test_that("PomaCreateObject handles special characters and whitespace in column names", {
  data(iris)
  metadata <- data.frame(sample_id = paste0("sample_", 1:150), group = iris$Species)
  features <- iris[, 1:4]
  colnames(features) <- c("ROCK1", "Rock2", "aKG@", "Petal.Width")
  se_object <- PomaCreateObject(metadata = metadata, features = features)
  expect_true(all(grepl("^[A-Za-z0-9_\\.]+$", rownames(SummarizedExperiment::assay(se_object)))))
})

