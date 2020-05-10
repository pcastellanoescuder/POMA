context("PomaClust")

test_that("PomaClust works", {
  
  data("st000284")
  data("st000336")
  imp_st000336 <- PomaImpute(st000336, method = "knn")
  
  a <- PomaClust(st000284)
  b <- PomaClust(imp_st000336)
  
  c <- PomaClust(st000284, method = "maximum", k = 5, show_clusters = FALSE, show_labels = TRUE)
  d <- PomaClust(imp_st000336, method = "manhattan", k = 2, show_clusters = FALSE, show_labels = TRUE)
  
  ## table
  
  expect_equal(nrow(a$mds_values), nrow(c$mds_values))
  expect_equal(nrow(b$mds_values), nrow(d$mds_values))
  
  expect_equal(4, ncol(a$mds_values))
  expect_equal(4, ncol(b$mds_values))
  expect_equal(4, ncol(c$mds_values))
  expect_equal(4, ncol(d$mds_values))
  
  ## plot
  
  expect_equal(class(a$mds_plot)[2], "ggplot")
  expect_equal(class(b$mds_plot)[2], "ggplot")
  expect_equal(class(c$mds_plot)[2], "ggplot")
  expect_equal(class(d$mds_plot)[2], "ggplot")
  
  ## errors
  
  expect_error(PomaClust())
  expect_error(PomaClust(iris))
  expect_error(PomaClust(st000284, method =  "euclid"))
  expect_error(PomaClust(st000284, method =  "max"))
  
})

