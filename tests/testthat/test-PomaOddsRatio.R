context("PomaOddsRatio")

test_that("PomaOddsRatio works", {

  data("st000284")

  norm_none <- PomaNorm(st000284, method = "none")
  norm_ls <- PomaNorm(st000284, method = "log_scaling")

  a <- PomaOddsRatio(norm_none)$OddsRatioPlot
  b <- PomaOddsRatio(norm_ls)$OddsRatioPlot

  c <- PomaOddsRatio(norm_ls, feature_name = "methyl_succinate_131_0_113_0")$OddsRatioPlot
  d <- PomaOddsRatio(norm_ls, feature_name = c("methyl_succinate_131_0_113_0", "linoleic_acid_277_1_259_0"))$OddsRatioPlot

  df_a <- layer_data(a)
  df_b <- layer_data(b)

  df_c <- layer_data(c)
  df_d <- layer_data(d)

  ##

  e <- PomaOddsRatio(norm_none)$OddsRatioTable
  f <- PomaOddsRatio(norm_ls)$OddsRatioTable
  g <- PomaOddsRatio(norm_ls, feature_name = "methyl_succinate_131_0_113_0")$OddsRatioTable
  h <- PomaOddsRatio(norm_ls, feature_name = c("methyl_succinate_131_0_113_0", "linoleic_acid_277_1_259_0"))$OddsRatioTable

  ##

  expect_equal(dim(e), dim(f))
  expect_equal(dim(f), dim(g))
  expect_equal(dim(g), dim(h))

  ##

  expect_equal(nrow(df_a), nrow(df_b))
  expect_false(nrow(df_c) == nrow(df_d))

  ##

  expect_error(PomaOddsRatio(norm_ls, feature_name = "hello"))
  expect_error(PomaOddsRatio(norm_ls, feature_name = "methyl_succinate_131_0_113_"))

})

