context("PomaBoxplots")

test_that("PomaBoxplots works", {
  
  data("st000284")
  
  norm_none <- PomaNorm(st000284, method = "none")
  norm_ls <- PomaNorm(st000284, method = "log_scaling")
  
  a <- PomaBoxplots(norm_none, label_size = 12)
  b <- PomaBoxplots(norm_ls, label_size = 10)
  c <- PomaBoxplots(norm_none, x = "features")
  d <- PomaBoxplots(norm_ls, x = "features")
  
  e <- PomaBoxplots(norm_none, x = "samples")
  
  f <- PomaBoxplots(norm_none, x = "samples", violin = TRUE)
  g <- PomaBoxplots(norm_none, x = "samples", violin = FALSE)
  h <- PomaBoxplots(norm_none, x = "features", violin = TRUE)
  i <- PomaBoxplots(norm_none, x = "features", violin = FALSE)
  
  j <- PomaBoxplots(norm_ls, x = "features", feature_name = "methyl_succinate")
  k <- PomaBoxplots(norm_ls, x = "features", feature_name = c("methyl_succinate", "linoleic_acid"))
  
  
  df_a <- ggplot2::layer_data(a)
  df_b <- ggplot2::layer_data(b)
  df_c <- ggplot2::layer_data(c)
  df_d <- ggplot2::layer_data(d)
  df_e <- ggplot2::layer_data(e)
  
  df_f <- ggplot2::layer_data(f)
  df_g <- ggplot2::layer_data(g)
  df_h <- ggplot2::layer_data(h)
  df_i <- ggplot2::layer_data(i)
  
  df_j <- ggplot2::layer_data(j)
  df_k <- ggplot2::layer_data(k)
  
  ####
  
  expect_true(min(df_a$ymin) == min(df_c$ymin))
  expect_true(min(df_b$ymin) != min(df_d$ymin))
  expect_false(all(df_a$ymin == df_b$ymin))
  expect_false(all(df_c$ymin == df_d$ymin))
  
  expect_equal(df_a, df_e)
  
  expect_false(nrow(df_j) == nrow(df_k))
  expect_false(nrow(df_j) == nrow(df_h))
  expect_false(nrow(df_k) == nrow(df_i))
  
  ##
  
  expect_error(PomaBoxplots(norm_ls, x = "samp"))
  expect_error(PomaBoxplots(x = "sample"))
  expect_error(PomaBoxplots(iris, x = "sample"))
  
  ##
  
  expect_error(PomaBoxplots(norm_ls, x = "features", feature_name = "hello"))
  expect_error(PomaBoxplots(norm_ls, x = "features", feature_name = "methyl_succina"))
  
  expect_message(PomaBoxplots(norm_ls, x = "features", feature_name = c("methyl_succina", "linoleic_acid")))
  expect_message(PomaBoxplots(norm_ls, x = "features", feature_name = c("methyl_succinate", "linoleic_aci")))
  
})

