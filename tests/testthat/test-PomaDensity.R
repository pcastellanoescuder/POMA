context("PomaDensity")

test_that("PomaDensity works", {

  data("st000284")

  norm_none <- PomaNorm(st000284, method = "none")
  norm_ls <- PomaNorm(st000284, method = "log_scaling")

  a <- PomaDensity(norm_none)
  b <- PomaDensity(norm_ls)
  c <- PomaDensity(norm_none, group = "features")
  d <- PomaDensity(norm_ls, group = "features")

  e <- PomaDensity(norm_none, group = "samples")

  f <- PomaDensity(norm_ls, group = "features", feature_name = "methyl_succinate")
  g <- PomaDensity(norm_ls, group = "features", feature_name = c("methyl_succinate", "linoleic_acid"))

  df_a <- ggplot2::layer_data(a)
  df_b <- ggplot2::layer_data(b)
  df_c <- ggplot2::layer_data(c)
  df_d <- ggplot2::layer_data(d)
  df_e <- ggplot2::layer_data(e)
  df_f <- ggplot2::layer_data(f)
  df_g <- ggplot2::layer_data(g)

  ####

  expect_false(all(df_a$y == df_c$y))
  expect_false(all(df_b$y == df_d$y))
  expect_false(all(df_a$y == df_b$y))
  expect_false(all(df_c$y == df_d$y))

  expect_false(all(df_f$y == df_g$y))
  expect_false(all(df_f$y == df_d$y))
  expect_false(all(df_g$y == df_d$y))

  expect_equal(df_a, df_e)

  expect_error(PomaDensity(norm_ls, group = "samp"))
  expect_error(PomaDensity(norm_ls, group = "features", feature_name = "hello"))

  expect_error(PomaDensity(norm_ls, feature_name = "hello"))
  expect_error(PomaDensity(norm_ls, feature_name = "methyl_succinat"))
  
  ##
  
  expect_error(PomaDensity(group = "sample"))
  expect_error(PomaDensity(iris, group = "sample"))
  expect_error(PomaDensity(norm_ls, group = "features", feature_name = c("methyl_succinate", "linoleic_acid"), 
                           legend_position = "no"))
  
  expect_warning(PomaDensity(norm_ls, group = "features", feature_name = c("methyl_succinate", "linoleic_aci")))
  expect_warning(PomaDensity(norm_ls, feature_name = c("methyl_succinat", "linoleic_acid")))
  
})

