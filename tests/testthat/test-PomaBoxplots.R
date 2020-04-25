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
  
  j <- PomaBoxplots(norm_ls, group = "features", feature_name = "methyl_succinate_131_0_113_0")
  k <- PomaBoxplots(norm_ls, group = "features", feature_name = c("methyl_succinate_131_0_113_0", "linoleic_acid_277_1_259_0"))
  
  
  df_a <- layer_data(a)
  df_b <- layer_data(b)
  df_c <- layer_data(c)
  df_d <- layer_data(d)
  df_e <- layer_data(e)
  
  df_f <- layer_data(f)
  df_g <- layer_data(g)
  df_h <- layer_data(h)
  df_i <- layer_data(i)
  
  df_j <- layer_data(j)
  df_k <- layer_data(k)
  
  ####
  
  expect_false(all(df_a$ymin == df_c$ymin))
  expect_false(all(df_b$ymin == df_d$ymin))
  expect_false(all(df_a$ymin == df_b$ymin))
  expect_false(all(df_c$ymin == df_d$ymin))
  
  expect_equal(df_f$outliers, df_g$outliers)
  expect_equal(df_h$outliers, df_i$outliers)
  
  expect_equal(df_a, df_e)
  
  expect_false(nrow(df_j) == nrow(df_k))
  expect_false(nrow(df_j) == nrow(df_h))
  expect_false(nrow(df_k) == nrow(df_i))
  
  ##
  
  expect_warning(PomaBoxplots(norm_none))
  expect_warning(PomaBoxplots(norm_ls))
  
  expect_error(PomaBoxplots(norm_ls, group = "samp"))
  
  ##
  
  expect_error(PomaBoxplots(group = "sample"))
  expect_error(PomaBoxplots(iris, group = "sample"))
  
  ##
  
  expect_error(PomaBoxplots(norm_ls, group = "features", feature_name = "hello"))
  expect_error(PomaBoxplots(norm_ls, group = "features", feature_name = "methyl_succinate_131_0_113_"))
  expect_error(PomaBoxplots(norm_ls, group = "features", feature_name = c("methyl_succinate_131_0_113_", "linoleic_acid_277_1_259_0")))
  
})

