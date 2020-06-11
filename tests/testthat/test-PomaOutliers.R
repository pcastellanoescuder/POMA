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
  
  ## TABLES
  
  expect_equal(nrow(a$outliers), nrow(Biobase::pData(st000284)) - nrow(Biobase::pData(b)))
  expect_equal(nrow(c$outliers), nrow(Biobase::pData(st000284)) - nrow(Biobase::pData(d)))
  
  expect_equal(nrow(e$outliers), nrow(Biobase::pData(st000336)) - nrow(Biobase::pData(f)))
  expect_equal(nrow(g$outliers), nrow(Biobase::pData(st000336)) - nrow(Biobase::pData(h)))
  
  expect_false(nrow(a$outliers) == nrow(c$outliers))
  expect_false(nrow(e$outliers) == nrow(g$outliers))
  
  ## PLOTS
  
  # df_a <- layer_data(lasso_res$coefficientPlot)
  # df_b <- layer_data(ridge_res$coefficientPlot)
  # df_e <- layer_data(enet_res$coefficientPlot)
  # df_c <- layer_data(lasso_res$cvLassoPlot)
  # 
  # expect_false(length(df_a$y) == length(df_b$y))
  # expect_false(length(df_a$y) == length(df_e$y))
  
  ## ERRORS
  
  # expect_error(PomaOutliers(alpha = 2))
  # expect_error(PomaOutliers(alpha = -0.5))
  # expect_error(PomaOutliers(iris))
  # expect_error(PomaOutliers(normalized_test, alpha = 1))
  # expect_error(PomaOutliers(normalized_test_less, alpha = 1))
  # expect_error(PomaOutliers())
  # expect_error(PomaOutliers(ntest = 60))
  # expect_error(PomaOutliers(ntest = 2))
  
})

