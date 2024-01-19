
test_that("PomaBatch handles valid SummarizedExperiment objects", {
  data <- create_mock_summarized_experiment()
  SummarizedExperiment::colData(data)$batch <- factor(rep(c("Batch1", "Batch2"), each = ncol(data)/2))
  corrected_data <- PomaBatch(data, batch = "batch")
  expect_is(corrected_data, "SummarizedExperiment")
})

test_that("PomaBatch stops with non-SummarizedExperiment objects", {
  data <- data.frame(matrix(runif(100), ncol = 10))
  expect_error(PomaBatch(data, batch = "batch"))
})

test_that("PomaBatch stops if batch column is not found", {
  data <- create_mock_summarized_experiment()
  expect_error(PomaBatch(data, batch = "nonexistent_batch"))
})

test_that("PomaBatch handles additional covariates", {
  data <- create_mock_summarized_experiment()
  SummarizedExperiment::colData(data)$batch <- factor(rep(c("Batch1", "Batch2"), each = ncol(data)/2))
  SummarizedExperiment::colData(data)$covariate <- rnorm(nrow(SummarizedExperiment::colData(data)))
  corrected_data <- PomaBatch(data, batch = "batch", mod = c("covariate"))
  expect_is(corrected_data, "SummarizedExperiment")
})

