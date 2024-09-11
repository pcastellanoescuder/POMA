
test_that("PomaRankProd handles valid SummarizedExperiment objects", {
  data <- create_mock_summarized_experiment(binary = TRUE)
  rank_prod_results <- PomaRankProd(data)
  expect_is(rank_prod_results, "list")
  expect_true(all(c("up_regulated", "down_regulated") %in% names(rank_prod_results)))
})

test_that("PomaRankProd stops with non-SummarizedExperiment objects", {
  data <- data.frame(matrix(runif(100), ncol = 10))
  expect_error(PomaRankProd(data), "data is not a SummarizedExperiment object")
})

test_that("PomaRankProd handles different parameters correctly", {
  data <- create_mock_summarized_experiment(binary = TRUE)
  rank_prod_results_logged <- PomaRankProd(data, logged = TRUE)
  rank_prod_results_cutoff <- PomaRankProd(data, cutoff = 0.1)
  rank_prod_results_method <- PomaRankProd(data, method = "pval")
  expect_is(rank_prod_results_logged, "list")
  expect_is(rank_prod_results_cutoff, "list")
  expect_is(rank_prod_results_method, "list")
})

test_that("PomaRankProd stops with incorrect method argument", {
  data <- create_mock_summarized_experiment(binary = TRUE)
  expect_error(PomaRankProd(data, method = "invalid_method"), "Incorrect value for method argument")
})

