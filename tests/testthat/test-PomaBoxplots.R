context("PomaBoxplots")

test_that("PomaBoxplots works", {
  
  # library(tidyverse)
  
  data("st000284")
  
  norm_none <- PomaNorm(st000284, method = "none")
  norm_ls <- PomaNorm(st000284, method = "log_scaling")
  
  a <- PomaBoxplots(norm_none)
  b <- PomaBoxplots(norm_ls)
  c <- PomaBoxplots(norm_none, group = "features")
  d <- PomaBoxplots(norm_ls, group = "features")
  
  e <- PomaBoxplots(norm_none, group = "samples")
  
  f <- PomaBoxplots(norm_none, group = "samples", jitter = T)
  g <- PomaBoxplots(norm_none, group = "samples", jitter = F)
  h <- PomaBoxplots(norm_none, group = "features", jitter = T)
  i <- PomaBoxplots(norm_none, group = "features", jitter = F)
  
  df_a <- layer_data(a)
  df_b <- layer_data(b)
  df_c <- layer_data(c)
  df_d <- layer_data(d)
  df_e <- layer_data(e)
  
  df_f <- layer_data(f)
  df_g <- layer_data(g)
  df_h <- layer_data(h)
  df_i <- layer_data(i)
  
  ####
  
  expect_false(all(df_a$ymin == df_c$ymin))
  expect_false(all(df_b$ymin == df_d$ymin))
  expect_false(all(df_a$ymin == df_b$ymin))
  expect_false(all(df_c$ymin == df_d$ymin))
  
  expect_equal(df_f$outliers, df_g$outliers)
  expect_equal(df_h$outliers, df_i$outliers)
  
  expect_equal(df_a, df_e)
  
  expect_warning(PomaBoxplots(norm_none))
  expect_warning(PomaBoxplots(norm_ls))
  
  expect_error(PomaBoxplots(norm_ls, group = "samp"))
  
})

