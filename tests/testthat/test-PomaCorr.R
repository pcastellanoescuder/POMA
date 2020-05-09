context("PomaCorr")

test_that("PomaCorr works", {
  
  data("st000284")
  data("st000336")
  imp_st000336 <- PomaImpute(st000336, method = "knn")
  
  a <- PomaCorr(st000284)
  b <- PomaCorr(imp_st000336)
  
  c <- PomaCorr(st000284, corr_type = "glasso", threshold = 0.3, rho = 0.9)
  d <- PomaCorr(imp_st000336 , corr_type = "glasso", threshold = 0.5, rho = 0.7)
  
  ## table
  
  expect_equal(((113*113)-113)/2, nrow(a$correlations))
  expect_equal(((30*30)-30)/2, nrow(b$correlations))
  expect_equal(((113*113)-113)/2, nrow(c$correlations))
  expect_equal(((30*30)-30)/2, nrow(d$correlations))
  
  expect_equal(3, ncol(a$correlations))
  expect_equal(3, ncol(b$correlations))
  expect_equal(3, ncol(c$correlations))
  expect_equal(3, ncol(d$correlations))
  
  ## corrplot
  
  expect_equal(class(a$corrplot)[2], "ggplot")
  expect_equal(class(b$corrplot)[2], "ggplot")
  expect_equal(class(c$corrplot)[2], "ggplot")
  expect_equal(class(d$corrplot)[2], "ggplot")
  
  ## networks
  
  expect_true(class(a$graph) == "qgraph")
  expect_true(class(b$graph) == "qgraph")
  expect_true(class(c$graph) == "qgraph")
  expect_true(class(d$graph) == "qgraph")
  
  expect_equal(113, nrow(a$graph$layout))
  expect_equal(30, nrow(b$graph$layout))
  expect_equal(113, nrow(c$graph$layout))
  expect_equal(30, nrow(d$graph$layout))
  
  ## errors
  
  expect_error(PomaCorr())
  expect_error(PomaCorr(iris))
  expect_error(PomaCorr(st000284, shape = "cir"))
  expect_error(PomaCorr(st000284, type = "lo"))
  expect_error(PomaCorr(st000284, corr_type = "co"))
  
})

