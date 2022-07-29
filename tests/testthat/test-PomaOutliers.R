context("PomaOutliers")

test_that("PomaOutliers works", {
  
  data("st000336")
  data("st000284")
  
  a <- st000284 %>% PomaImpute() %>% PomaNorm() %>% PomaOutliers(do = "analyze")
  b <- st000284 %>% PomaImpute() %>% PomaNorm() %>% PomaOutliers()
  c <- st000284 %>% PomaImpute() %>% PomaNorm() %>% PomaOutliers(do = "analyze", coef = 3)
  d <- st000284 %>% PomaImpute() %>% PomaNorm() %>% PomaOutliers(coef = 3)

  e <- st000336 %>% PomaImpute() %>% PomaNorm() %>% PomaOutliers(do = "analyze")
  f <- st000336 %>% PomaImpute() %>% PomaNorm() %>% PomaOutliers()
  g <- st000336 %>% PomaImpute() %>% PomaNorm() %>% PomaOutliers(do = "analyze", coef = 3)
  h <- st000336 %>% PomaImpute() %>% PomaNorm() %>% PomaOutliers(coef = 3)
  
  i <- st000336 %>% PomaImpute() %>% PomaNorm() %>% PomaOutliers(do = "analyze", method = "maximum", type = "centroid", labels = TRUE)
  j <- st000336 %>% PomaImpute() %>% PomaNorm() %>% PomaOutliers(do = "analyze", method = "manhattan", type = "centroid", labels = FALSE)
  k <- st000336 %>% PomaImpute() %>% PomaNorm() %>% PomaOutliers(do = "analyze", method = "canberra", type = "centroid", labels = TRUE)
  l <- st000336 %>% PomaImpute() %>% PomaNorm() %>% PomaOutliers(do = "analyze", method = "minkowski", type = "centroid", labels = FALSE)

  ##
  
  expect_equal(nrow(a$outliers), nrow(SummarizedExperiment::colData(st000284)) - nrow(SummarizedExperiment::colData(b)))
  expect_equal(nrow(c$outliers), nrow(SummarizedExperiment::colData(st000284)) - nrow(SummarizedExperiment::colData(d)))
  
  expect_equal(nrow(e$outliers), nrow(SummarizedExperiment::colData(st000336)) - nrow(SummarizedExperiment::colData(f)))
  expect_equal(nrow(g$outliers), nrow(SummarizedExperiment::colData(st000336)) - nrow(SummarizedExperiment::colData(h)))
  
  expect_false(nrow(a$outliers) == nrow(c$outliers))
  expect_false(nrow(e$outliers) == nrow(g$outliers))
  
  expect_equal(ncol(i$outliers), ncol(j$outliers))
  expect_equal(ncol(j$outliers), ncol(k$outliers))
  expect_equal(ncol(k$outliers), ncol(l$outliers))
  
  ## PLOTS
  
  df_i <- ggplot2::layer_data(i$polygon_plot)
  df_j <- ggplot2::layer_data(j$polygon_plot)
  
  df_k <- ggplot2::layer_data(k$distance_boxplot)
  df_l <- ggplot2::layer_data(l$distance_boxplot)
  
  expect_false(length(df_i$y) == length(df_j$y))
  expect_equal(ncol(df_i), ncol(df_j))
  
  expect_false(all(df_k$ymin == df_l$ymin))
  expect_equal(ncol(df_k), ncol(df_l))
  
  ## ERRORS
  
  expect_error(PomaOutliers())
  expect_error(PomaOutliers(iris))
  expect_error(PomaOutliers(st000284, method = "man"))
  expect_error(PomaOutliers(st000284, type = "media"))
  expect_error(PomaOutliers(st000284, do = "analise"))
  
})

