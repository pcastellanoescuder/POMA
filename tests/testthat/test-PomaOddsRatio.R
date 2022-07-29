context("PomaOddsRatio")

test_that("PomaOddsRatio works", {

  data("st000336")

  imputed <- POMA::PomaImpute(st000336, method = "knn")
  
  norm_none <- PomaNorm(imputed, method = "none")
  norm_ls <- PomaNorm(imputed, method = "log_scaling")

  a <- PomaOddsRatio(norm_none, feature_name = c("glutamic_acid", "glutamine", "glycine", "histidine", "isoleucine", "leucine", "lysine"))$OddsRatioPlot
  b <- PomaOddsRatio(norm_ls, feature_name = c("glutamic_acid", "glutamine", "glycine", "histidine", "isoleucine", "leucine", "lysine"))$OddsRatioPlot

  c <- PomaOddsRatio(norm_ls, feature_name = "glutamic_acid")$OddsRatioPlot
  d <- PomaOddsRatio(norm_ls, feature_name = c("glutamic_acid", "arginine"))$OddsRatioPlot

  df_a <- ggplot2::layer_data(a)
  df_b <- ggplot2::layer_data(b)

  ##

  expect_equal(nrow(df_a), nrow(df_b))

  ##
  
  e <- PomaOddsRatio(norm_none)$OddsRatioTable
  f <- PomaOddsRatio(norm_ls)$OddsRatioTable
  g <- PomaOddsRatio(norm_ls, feature_name = "glutamic_acid")$OddsRatioTable
  h <- PomaOddsRatio(norm_ls, feature_name = c("glutamic_acid", "arginine"))$OddsRatioTable
  h_1 <- PomaOddsRatio(norm_ls, feature_name = c("glutamic_acid", "arginine"), covariates = TRUE)$OddsRatioTable

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

  k <- PomaOddsRatio(norm_ls, feature_name = c("glutamic_acid", "glutamine", "glycine", "histidine", "isoleucine", "leucine", "lysine"))$OddsRatioPlot
  l <- PomaOddsRatio(norm_ls, feature_name = c("glutamic_acid", "glutamine", "glycine", "histidine", "isoleucine", "leucine", "lysine"), showCI = FALSE)$OddsRatioPlot
  m <- PomaOddsRatio(norm_ls, feature_name = c("glutamic_acid", "glutamine", "glycine", "histidine", "isoleucine", "leucine", "lysine"), showCI = FALSE, covariates = TRUE)$OddsRatioPlot

  df_k <- ggplot2::layer_data(k)
  df_l <- ggplot2::layer_data(l)
  df_m <- ggplot2::layer_data(m)

  ##

  expect_equal(nrow(df_k), nrow(df_l))

  ##

  expect_error(PomaOddsRatio(norm_ls, feature_name = "hello"))
  expect_error(PomaOddsRatio(norm_ls, feature_name = "glutamic_aci"))
  expect_error(PomaOddsRatio(norm_ls, feature_name = c("glutamic_aci", "arginine")))

  ##

  SummarizedExperiment::colData(imputed) <- SummarizedExperiment::colData(imputed)[1]
  expect_error(PomaOddsRatio(imputed, covariates = TRUE))
  
  ##
  
  expect_error(PomaOddsRatio())
  expect_error(PomaOddsRatio(iris))

})

