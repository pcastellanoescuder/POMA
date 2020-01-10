context("PomaRandForest")

test_that("PomaRandForest works", {

  # library(tidyverse)

  data("st000284")

  res <- PomaRandForest(st000284)
  res1 <- PomaRandForest(st000284)

  df_a <- layer_data(res$error_tree)
  df_b <- layer_data(res$gini_plot)

  df_c <- layer_data(res1$error_tree)
  df_d <- layer_data(res1$gini_plot)

  ####

  expect_error(PomaRandForest())
  expect_true(length(res) == 5)

  expect_false(all(df_a$y == df_c$y))
  expect_false(all(df_b$y == df_d$y))

  expect_equal(length(df_a$y), length(df_c$y))
  expect_equal(length(df_b$y), length(df_d$y))

  expect_false(all(res$importance_pred == res1$importance_pred))
  expect_equal(nrow(res$importance_pred), nrow(res1$importance_pred))

  expect_false(all(res$forest_data == res1$forest_data))
  expect_equal(nrow(res$forest_data), nrow(res1$forest_data))

  expect_equal(dim(res$confusion_matrix), dim(res1$confusion_matrix))
  expect_false(all(res$confusion_matrix == res1$confusion_matrix))

})

