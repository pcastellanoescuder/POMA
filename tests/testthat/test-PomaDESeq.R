
test_that("PomaDESeq handles valid SummarizedExperiment objects", {
  data <- create_mock_summarized_experiment(integers = TRUE)
  DESeq_results <- PomaDESeq(data)
  expect_is(DESeq_results, "tbl_df")
})

test_that("PomaDESeq stops with non-SummarizedExperiment objects", {
  data <- data.frame(matrix(runif(100), ncol = 10))
  expect_error(PomaDESeq(data), "data is not a SummarizedExperiment object")
})

test_that("PomaDESeq handles different adjust methods correctly", {
  data <- create_mock_summarized_experiment(integers = TRUE)
  for (adjust_method in c("fdr", "holm", "hochberg", "hommel", "bonferroni", "BH", "BY")) {
    DESeq_results <- PomaDESeq(data, adjust = adjust_method)
    expect_is(DESeq_results, "tbl_df")
  }
})

test_that("PomaDESeq requires metadata", {
  data <- create_mock_summarized_experiment(integers = TRUE)
  metadata_removed_data <- SummarizedExperiment::SummarizedExperiment(assays = SummarizedExperiment::assay(data))
  expect_error(PomaDESeq(metadata_removed_data), "metadata file required")
})

