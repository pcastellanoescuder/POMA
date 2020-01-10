context("PomaVolcano")

test_that("PomaVolcano works", {

  library(tidyverse)
  library(ggrepel)

  data("st000284")

  a <- PomaVolcano(st000284)
  b <- PomaVolcano(st000284, pval = "adjusted", pval_cutoff = 0.05, log2FC = 0.6, xlim = 2)
  c <- PomaVolcano(st000284, pval = "raw", pval_cutoff = 0.05, log2FC = 0.6, xlim = 2)

  df_a <- layer_data(a)
  df_b <- layer_data(b)
  df_c <- layer_data(c)

  expect_false(all(df_a$y == df_b$y))
  expect_equal(df_a, df_c)

  expect_warning(PomaVolcano(st000284))
  expect_error(PomaVolcano(st000284, pval = "ra"))

})
