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

  f <- PomaDensity(norm_ls, group = "features", feature_name = "methyl_succinate_131_0_113_0")
  g <- PomaDensity(norm_ls, group = "features", feature_name = c("methyl_succinate_131_0_113_0", "linoleic_acid_277_1_259_0"))

  df_a <- layer_data(a)
  df_b <- layer_data(b)
  df_c <- layer_data(c)
  df_d <- layer_data(d)
  df_e <- layer_data(e)
  df_f <- layer_data(f)
  df_g <- layer_data(g)

  ####

  expect_false(all(df_a$y == df_c$y))
  expect_false(all(df_b$y == df_d$y))
  expect_false(all(df_a$y == df_b$y))
  expect_false(all(df_c$y == df_d$y))

  expect_false(all(df_f$y == df_g$y))
  expect_false(all(df_f$y == df_d$y))
  expect_false(all(df_g$y == df_d$y))

  expect_equal(df_a, df_e)

  expect_warning(PomaDensity(norm_none))
  expect_warning(PomaDensity(norm_ls))

  expect_error(PomaDensity(norm_ls, group = "samp"))
  expect_error(PomaDensity(norm_ls, group = "features", feature_name = "hello"))

  expect_error(PomaDensity(norm_ls, feature_name = "hello"))
  expect_error(PomaDensity(norm_ls, feature_name = "methyl_succinate_131_0_113_"))
  expect_error(PomaDensity(norm_ls, feature_name = c("methyl_succinate_131_0_113_", "linoleic_acid_277_1_259_0")))
  
  ##
  
  expect_error(PomaDensity(group = "sample"))
  expect_error(PomaDensity(iris, group = "sample"))
  
})

