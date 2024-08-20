
test_that("PomaVolcano handles valid data frames", {
  data <- create_mock_summarized_experiment(binary = TRUE) %>% 
    PomaImpute() %>% 
    PomaUnivariate() %>%
    dplyr::select(feature, fold_change, pvalue)
  plot <- PomaVolcano(data)
  expect_is(plot, "ggplot")
})

test_that("PomaVolcano handles different pval_cutoff and log2fc_cutoff values", {
  data <- create_mock_summarized_experiment(binary = TRUE) %>% 
    PomaImpute() %>% 
    PomaUnivariate() %>%
    dplyr::select(feature, fold_change, pvalue)
  plot <- PomaVolcano(data, pval_cutoff = 0.01, log2fc_cutoff = 1)
  expect_is(plot, "ggplot")
})

test_that("PomaVolcano handles labels parameter correctly", {
  data <- create_mock_summarized_experiment(binary = TRUE) %>% 
    PomaImpute() %>% 
    PomaUnivariate() %>%
    dplyr::select(feature, fold_change, pvalue)
  plot_with_labels <- PomaVolcano(data, labels = TRUE)
  expect_is(plot_with_labels, "ggplot")
})

test_that("PomaVolcano handles custom axis labels", {
  data <- create_mock_summarized_experiment(binary = TRUE) %>% 
    PomaImpute() %>% 
    PomaUnivariate() %>%
    dplyr::select(feature, fold_change, pvalue)
  plot_custom_labels <- PomaVolcano(data, x_label = "Custom X", y_label = "Custom Y")
  expect_is(plot_custom_labels, "ggplot")
})

