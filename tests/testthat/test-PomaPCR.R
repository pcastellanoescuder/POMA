
test_that("PomaPCR handles valid SummarizedExperiment objects", {
  data <- create_mock_summarized_experiment()
  data$NumericVariable <- 1:20
  pcr_results <- PomaPCR(data, y = "NumericVariable")
  expect_is(pcr_results, "tbl_df")
  expect_true(all(c("component", "estimate", "std_err", "statistic", "pvalue", "adj_pvalue") %in% names(pcr_results)))
})

test_that("PomaPCR stops with non-SummarizedExperiment objects", {
  data <- data.frame(matrix(runif(100), ncol = 10))
  expect_error(PomaPCR(data, y = "NumericVariable"), "data is not a SummarizedExperiment object")
})

test_that("PomaPCR handles different number of components", {
  data <- create_mock_summarized_experiment()
  data$NumericVariable <- 1:20
  pcr_results_2comp <- PomaPCR(data, ncomp = 2, y = "NumericVariable")
  pcr_results_5comp <- PomaPCR(data, ncomp = 5, y = "NumericVariable")
  expect_is(pcr_results_2comp, "tbl_df")
  expect_is(pcr_results_5comp, "tbl_df")
})

test_that("PomaPCR stops with incorrect y variable", {
  data <- create_mock_summarized_experiment()
  data$InvalidVariable <- rep("a", 20)
  expect_error(PomaPCR(data, y = "InvalidVariable"), "No numeric variables to be used as dependent variable in metadata file")
})

test_that("PomaPCR stops with incorrect adjust argument", {
  data <- create_mock_summarized_experiment()
  data$NumericVariable <- 1:20
  expect_error(PomaPCR(data, y = "NumericVariable", adjust = "invalid"), "Incorrect value for adjust argument")
})

