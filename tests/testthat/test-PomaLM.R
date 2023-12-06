
test_that("PomaLM handles valid SummarizedExperiment objects", {
  data <- create_mock_summarized_experiment()
  data$NumericVariable <- 1:20
  lm_results <- PomaLM(data)
  expect_is(lm_results, "list")
  expect_true("lm_table" %in% names(lm_results))
  expect_true("regression_plot" %in% names(lm_results))
})

test_that("PomaLM stops with non-SummarizedExperiment objects", {
  data <- data.frame(matrix(runif(100), ncol = 10))
  expect_error(PomaLM(data), "data is not a SummarizedExperiment object")
})

test_that("PomaLM handles specific independent variables", {
  data <- create_mock_summarized_experiment()
  data$NumericVariable <- 1:20
  lm_results <- PomaLM(data, x = c("V1", "V2"))
  expect_is(lm_results, "list")
})

test_that("PomaLM handles specific dependent variables", {
  data <- create_mock_summarized_experiment()
  data$NumericVariable <- 1:20
  lm_results <- PomaLM(data, y = "NumericVariable")
  expect_is(lm_results, "list")
})

test_that("PomaLM provides expected output structure", {
  data <- create_mock_summarized_experiment()
  data$NumericVariable <- 1:20
  lm_results <- PomaLM(data)
  expect_is(lm_results$lm_table, "tbl_df")
  expect_true(all(c("feature", "estimate", "std_err", "statistic", "pvalue", "adj_pvalue") %in% names(lm_results$lm_table)))
})

