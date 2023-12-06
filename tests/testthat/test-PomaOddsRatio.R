
test_that("PomaOddsRatio handles valid SummarizedExperiment objects", {
  data <- create_mock_summarized_experiment(binary = TRUE)
  odds_ratio_results <- PomaOddsRatio(data)
  expect_is(odds_ratio_results, "list")
  expect_true(all(c("odds_ratio_table", "odds_ratio_plot") %in% names(odds_ratio_results)))
})

test_that("PomaOddsRatio stops with non-SummarizedExperiment objects", {
  data <- data.frame(matrix(runif(100), ncol = 10))
  expect_error(PomaOddsRatio(data), "data is not a SummarizedExperiment object")
})

test_that("PomaOddsRatio handles specific feature names correctly", {
  data <- create_mock_summarized_experiment(binary = TRUE)
  feature_name <- rownames(SummarizedExperiment::assay(data))[1:2]
  odds_ratio_results <- PomaOddsRatio(data, feature_name = feature_name)
  expect_is(odds_ratio_results, "list")
  expect_true(all(feature_name %in% odds_ratio_results$odds_ratio_table$feature))
})

test_that("PomaOddsRatio stops with incorrect feature names", {
  data <- create_mock_summarized_experiment(binary = TRUE)
  expect_error(PomaOddsRatio(data, feature_name = "non_existing_feature"))
})

test_that("PomaOddsRatio handles covariates and show_ci parameters correctly", {
  data <- create_mock_summarized_experiment(binary = TRUE)
  odds_ratio_results_with_covs <- PomaOddsRatio(data, covs = c("Group"))
  odds_ratio_results_no_ci <- PomaOddsRatio(data, show_ci = FALSE)
  expect_is(odds_ratio_results_with_covs, "list")
  expect_is(odds_ratio_results_no_ci, "list")
})

