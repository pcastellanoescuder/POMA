context("PomaHeatmap")

test_that("PomaHeatmap works", {
  
  data("st000284")
  data("st000336")
  
  a <- PomaHeatmap(st000284, sample_names = TRUE, feature_names = FALSE, show_legend = TRUE)
  b <- PomaHeatmap(st000284, sample_names = FALSE, feature_names = FALSE, show_legend = FALSE)
  
  c <- PomaHeatmap(st000336, sample_names = TRUE, feature_names = TRUE, show_legend = TRUE)
  d <- PomaHeatmap(st000336, sample_names = FALSE, feature_names = TRUE, show_legend = FALSE)
  
  ##
  
  expect_equal(class(a), class(b))
  expect_equal(class(b), class(d))
  
  expect_false(length(a@matrix) == length(d@matrix))
  
  expect_equal(length(a@row_order), ncol(t(MSnbase::exprs(st000284))))
  expect_equal(length(d@row_order), ncol(t(MSnbase::exprs(st000336))))
  
  ##
  
  expect_error(PomaHeatmap())
  expect_error(PomaHeatmap(iris))
  
})

