
test_that("PomaVolcano handles valid SummarizedExperiment objects", {
  data <- create_mock_summarized_experiment(binary = TRUE)
  plot <- PomaVolcano(data)
  expect_is(plot, "ggplot")
})

test_that("PomaVolcano stops with non-SummarizedExperiment objects", {
  data <- data.frame(matrix(runif(100), ncol = 10))
  expect_error(PomaVolcano(data), "data is not a SummarizedExperiment object")
})

test_that("PomaVolcano handles different methods correctly", {
  data <- create_mock_summarized_experiment(binary = TRUE)
  for (method in c("ttest", "mann", "limma")) { # DESeq
    plot <- PomaVolcano(data, method = method)
    expect_is(plot, "ggplot")
  }
})

test_that("PomaVolcano stops with incorrect method argument", {
  data <- create_mock_summarized_experiment(binary = TRUE)
  expect_error(PomaVolcano(data, method = "invalid_method"), "Incorrect value for method argument")
})

test_that("PomaVolcano handles p-value adjustments correctly", {
  data <- create_mock_summarized_experiment(binary = TRUE)
  for (adjust in c("fdr", "holm", "hochberg", "hommel", "bonferroni", "BH", "BY")) {
    plot <- PomaVolcano(data, adjust = adjust)
    expect_is(plot, "ggplot")
  }
})

test_that("PomaVolcano stops with incorrect pval argument", {
  data <- create_mock_summarized_experiment(binary = TRUE)
  expect_error(PomaVolcano(data, pval = "invalid_pval"), 'Incorrect value for pval argument')
})

test_that("PomaVolcano handles different pval_cutoff and log2fc_cutoff values", {
  data <- create_mock_summarized_experiment(binary = TRUE)
  plot <- PomaVolcano(data, pval_cutoff = 0.01, log2fc_cutoff = 1)
  expect_is(plot, "ggplot")
})

test_that("PomaVolcano handles labels, paired, and var_equal parameters correctly", {
  data <- create_mock_summarized_experiment(binary = TRUE)
  plot_with_labels <- PomaVolcano(data, labels = TRUE)
  plot_paired <- PomaVolcano(data, paired = TRUE)
  plot_var_equal <- PomaVolcano(data, var_equal = TRUE)
  expect_is(plot_with_labels, "ggplot")
  expect_is(plot_paired, "ggplot")
  expect_is(plot_var_equal, "ggplot")
})

