context("PomaHeatmap")

test_that("PomaHeatmap works", {
  
  data("st000284")
  data("st000336")
  
  a <- PomaHeatmap(st000284, scale = TRUE, scale_by = "features")
  b <- PomaHeatmap(st000284, scale = FALSE, scale_by = "features")
  
  c <- PomaHeatmap(st000336, scale = TRUE, scale_by = "samples")
  d <- PomaHeatmap(st000336, scale = FALSE, scale_by = "features")
  
  ##
  
  expect_equal(class(a), class(b))
  expect_equal(class(b), class(d))
  
  expect_false(length(a$rowInd) == length(d$rowInd))
  
  expect_equal(length(a$rowInd), ncol(t(Biobase::exprs(st000284))))
  expect_equal(length(d$rowInd), ncol(t(Biobase::exprs(st000336))))
  
  ##
  
  expect_error(PomaHeatmap())
  expect_error(PomaHeatmap(iris))
  expect_error(PomaHeatmap(st000284, scale_by = "samp"))
  
})

