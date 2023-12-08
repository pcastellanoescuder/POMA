
test_that("PomaLimma handles valid SummarizedExperiment objects", {
  data <- create_mock_summarized_experiment()
  limma_results <- PomaLimma(data, contrast = "A-B")
  expect_is(limma_results, "tbl_df")
})

test_that("PomaLimma stops with non-SummarizedExperiment objects", {
  data <- data.frame(matrix(runif(100), ncol = 10))
  expect_error(PomaLimma(data, contrast = "A-B"), "data is not a SummarizedExperiment object")
})

test_that("PomaLimma handles different contrasts correctly", {
  data <- create_mock_summarized_experiment()
  contrast <- levels(SummarizedExperiment::colData(data)[,1])
  limma_results <- PomaLimma(data, contrast = paste(contrast[1], contrast[2], sep = "-"))
  expect_is(limma_results, "tbl_df")
})

test_that("PomaLimma handles parameters correctly", {
  data <- create_mock_summarized_experiment()
  limma_results_with_covs <- PomaLimma(data, contrast = "A-B")
  limma_results_adjusted <- PomaLimma(data, contrast = "A-B", adjust = "holm")
  expect_is(limma_results_with_covs, "tbl_df")
  expect_is(limma_results_adjusted, "tbl_df")
})

