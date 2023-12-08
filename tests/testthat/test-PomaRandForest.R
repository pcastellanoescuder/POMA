
test_that("PomaRandForest handles valid SummarizedExperiment objects", {
  data <- create_mock_summarized_experiment()
  result <- PomaRandForest(data)
  expect_is(result, "list")
  expect_true(all(c("MeanDecreaseGini", "MeanDecreaseGini_plot", "oob_error", "error_tree", "model") %in% names(result)))
})

test_that("PomaRandForest stops with non-SummarizedExperiment objects", {
  data <- data.frame(matrix(runif(100), ncol = 10))
  expect_error(PomaRandForest(data), "data is not a SummarizedExperiment object")
})

test_that("PomaRandForest handles different ntest values correctly", {
  data <- create_mock_summarized_experiment()
  for (ntest in c(5, 10, 20, 30, 40, 50)) {
    result <- PomaRandForest(data, ntest = ntest)
    expect_is(result, "list")
    expect_true(all(c("confusionMatrix", "train_x", "train_y", "test_x", "test_y") %in% names(result)))
  }
})

test_that("PomaRandForest stops with invalid ntest argument", {
  data <- create_mock_summarized_experiment()
  expect_error(PomaRandForest(data, ntest = 60), "Incorrect value for ntest argument")
  expect_error(PomaRandForest(data, ntest = 0), "Incorrect value for ntest argument")
})

test_that("PomaRandForest handles different values for ntree, mtry, and nodesize", {
  data <- create_mock_summarized_experiment()
  result <- PomaRandForest(data, ntree = 100, mtry = 2, nodesize = 2)
  expect_is(result, "list")
})

test_that("PomaRandForest returns expected results with nvar parameter", {
  data <- create_mock_summarized_experiment()
  result <- PomaRandForest(data, nvar = 5)
  expect_is(result, "list")
  expect_equal(nrow(result$MeanDecreaseGini), 5)
})

