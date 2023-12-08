
test_that("PomaLMM handles valid SummarizedExperiment objects", {
  data <- create_mock_summarized_experiment()
  lmm_results <- PomaLMM(data)
  expect_is(lmm_results, "list")
  expect_true("variances" %in% names(lmm_results))
  expect_true("variances_plot" %in% names(lmm_results))
})

test_that("PomaLMM stops with non-SummarizedExperiment objects", {
  data <- data.frame(matrix(runif(100), ncol = 10))
  expect_error(PomaLMM(data), "data is not a SummarizedExperiment object")
})

test_that("PomaLMM handles specific independent variables", {
  data <- create_mock_summarized_experiment()
  data$NumericVariable <- 1:20
  data$FactorVariable <- factor(c(rep("s", 10), rep("d", 10)))
  lmm_results <- PomaLMM(data, x = c("NumericVariable", "FactorVariable"))
  expect_is(lmm_results, "list")
})

test_that("PomaLMM handles specific dependent variables", {
  data <- create_mock_summarized_experiment()
  lmm_results <- PomaLMM(data, y = c("V1", "V2"))
  expect_is(lmm_results, "list")
})

test_that("PomaLMM stops with incorrect adjustment method", {
  data <- create_mock_summarized_experiment()
  expect_error(PomaLMM(data, adjust = "invalid"), "Incorrect value for adjust argument")
})

test_that("PomaLMM provides expected output structure", {
  data <- create_mock_summarized_experiment()
  lmm_results <- PomaLMM(data)
  expect_is(lmm_results$variances, "tbl_df")
  expect_true(all(c("feature", "Residual", "(Intercept)") %in% names(lmm_results$variances)))
})

