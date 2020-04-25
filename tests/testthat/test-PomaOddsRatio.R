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

  expect_equal(nrow(df_a), nrow(df_b))
  expect_false(nrow(df_c) == nrow(df_d))

  ##

  e <- PomaOddsRatio(norm_none)$OddsRatioTable
  f <- PomaOddsRatio(norm_ls)$OddsRatioTable
  g <- PomaOddsRatio(norm_ls, feature_name = "methyl_succinate_131_0_113_0")$OddsRatioTable
  h <- PomaOddsRatio(norm_ls, feature_name = c("methyl_succinate_131_0_113_0", "linoleic_acid_277_1_259_0"))$OddsRatioTable
  h_1 <- PomaOddsRatio(norm_ls, feature_name = c("methyl_succinate_131_0_113_0", "linoleic_acid_277_1_259_0"), covariates = TRUE)$OddsRatioTable

  ##

  expect_equal(dim(e), dim(f))
  expect_false(nrow(f) == nrow(g))
  expect_false(nrow(g) == nrow(h))
  expect_false(nrow(h) == nrow(h_1))

  ##

  i <- PomaOddsRatio(norm_ls, feature_name = NULL, covariates = TRUE)$OddsRatioTable
  j <- PomaOddsRatio(norm_ls, feature_name = NULL, covariates = FALSE)$OddsRatioTable

  ##

  expect_false(nrow(i) == nrow(j))

  ##

  k <- PomaOddsRatio(norm_ls, feature_name = NULL)$OddsRatioPlot
  l <- PomaOddsRatio(norm_ls, feature_name = NULL, showCI = FALSE)$OddsRatioPlot
  m <- PomaOddsRatio(norm_ls, feature_name = NULL, showCI = FALSE, covariates = TRUE)$OddsRatioPlot

  df_k <- layer_data(k)
  df_l <- layer_data(l)
  df_m <- layer_data(m)

  ##

  expect_equal(nrow(df_k), nrow(df_l))
  expect_false(nrow(df_l) == nrow(df_m))

  ##

  expect_error(PomaOddsRatio(norm_ls, feature_name = "hello"))
  expect_error(PomaOddsRatio(norm_ls, feature_name = "methyl_succinate_131_0_113_"))
  expect_error(PomaOddsRatio(norm_ls, feature_name = c("methyl_succinate_131_0_113_", "linoleic_acid_277_1_259_0")))

  ##

  Biobase::pData(st000284) <- Biobase::pData(st000284)[1]
  expect_error(PomaOddsRatio(st000284, covariates = TRUE))
  
  ##
  
  expect_error(PomaOddsRatio())
  expect_error(PomaOddsRatio(iris))

})

