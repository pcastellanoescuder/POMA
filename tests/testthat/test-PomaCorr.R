context("PomaCorr")

test_that("PomaCorr works", {
  
  data("st000284")
  data("st000336")
  imp_st000336 <- PomaImpute(st000336, method = "knn")
  
  a <- PomaCorr(st000284)
  b <- PomaCorr(imp_st000336)
  
  c <- PomaCorr(st000284, corr_type = "glasso", coeff = 0.3, method = "spearman")
  d <- PomaCorr(imp_st000336 , corr_type = "glasso", coeff = 0.5, type = "upper")
  
  e <- PomaCorr(st000284, corr_type = "glasso", coeff = 0.5)
  f <- PomaCorr(st000284, coeff = 0.5)
  
  ## table
  
  expect_equal(((113*113)-113)/2, nrow(a$correlations))
  expect_equal(((30*30)-30)/2, nrow(b$correlations))
  expect_equal(((113*113)-113)/2, nrow(c$correlations))
  expect_equal(((30*30)-30)/2, nrow(d$correlations))
  
  expect_equal(3, ncol(a$correlations))
  expect_equal(3, ncol(b$correlations))
  expect_equal(3, ncol(c$correlations))
  expect_equal(3, ncol(d$correlations))
  
  expect_equal(class(e$data_glasso), "data.frame")
  expect_equal(ncol(e$data_glasso), 3)
  expect_equal(ncol(e$data_glasso), ncol(a$correlations))
  expect_equal(class(f$data_glasso), "NULL")
  
  ## corrplot
  
  expect_equal(class(a$corrplot)[2], "ggplot")
  expect_equal(class(b$corrplot)[2], "ggplot")
  expect_equal(class(c$corrplot)[2], "ggplot")
  expect_equal(class(d$corrplot)[2], "ggplot")
  
  ## networks
  
  expect_true(class(a$graph)[1] == "ggraph")
  expect_true(class(b$graph)[1] == "ggraph")
  expect_true(class(c$graph)[1] == "ggraph")
  expect_true(class(d$graph)[1] == "ggraph")
  
  expect_true(113 > nrow(a$graph$data))
  expect_equal(6, ncol(a$graph$data))
  
  expect_true(30 > nrow(b$graph$data))
  expect_equal(6, ncol(b$graph$data))
  
  expect_true(113 > nrow(c$graph$data))
  expect_true(30 > nrow(d$graph$data))
  
  ## errors
  
  expect_error(PomaCorr())
  expect_error(PomaCorr(iris))
  expect_error(PomaCorr(st000284, shape = "cir"))
  expect_error(PomaCorr(st000284, type = "lo"))
  expect_error(PomaCorr(st000284, corr_type = "co"))
  expect_error(PomaCorr(st000284, coeff = 2))
  expect_error(PomaCorr(st000284, coeff = -0.2))
  expect_error(PomaCorr(st000284, method = "pear"))
  expect_error(PomaCorr(st000284, corr_type = "cor", coeff = 1))
  expect_error(PomaCorr(st000284, corr_type = "glasso", coeff = 1))
  
})

