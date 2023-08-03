context("PomaRandForest")

test_that("PomaRandForest works", {

  # library(tidyverse)

  data("st000284")

  res <- PomaRandForest(st000284, ntest = 15, nvar = 15)
  res1 <- PomaRandForest(st000284, ntest = 30, mtry = 5)
  res2 <- PomaRandForest(st000284, mtry = 5)
  
  df_a <- ggplot2::layer_data(res$error_tree)
  df_b <- ggplot2::layer_data(res$MeanDecreaseGini_plot)

  df_c <- ggplot2::layer_data(res1$error_tree)
  df_d <- ggplot2::layer_data(res1$MeanDecreaseGini_plot)

  ####

  expect_error(PomaRandForest())
  expect_error(PomaRandForest(iris))
  expect_true(length(res) == 10)
  expect_false(length(res) == length(res2))
  expect_true(length(res2) == 5)
  
  expect_false(all(df_a$y == df_c$y))
  expect_false(length(df_b$y) == length(df_d$y))

  expect_equal(length(df_a$y), length(df_c$y))
  expect_false(length(df_b$y) == length(df_d$y))

  expect_false(all(res$oob_error == res1$oob_error))
  expect_equal(nrow(res$oob_error), nrow(res1$oob_error))

  expect_equal(dim(res$confusionMatrix), dim(res1$confusionMatrix))
  expect_false(all(res$confusionMatrix$table == res1$confusionMatrix$table))
  expect_equal(dim(res2$confusionMatrix), NULL)
  
  expect_error(PomaRandForest(st000284, ntest = 1))
  expect_error(PomaRandForest(st000284, ntest = 51))

})

