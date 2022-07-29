context("PomaImpute")

test_that("PomaImpute works", {
  
  library(SummarizedExperiment)
  
  data("st000284")

  data <- t(SummarizedExperiment::assay(st000284))

  data <- data*round(runif(n = 1, min = 0.01, max = 0.99), 3) # just to create decimals
  data[1:4, 5] <- 0 # create some zeros in the first group
  data[73:77, 5] <- 0 # create some zeros in the second group

  data[10:14, 5] <- NA # create some NA in the first group
  data[78:81, 5] <- NA # create some NA in the second group

  colnames(data) <- gsub("_", ":", colnames(data))
    
  target <- SummarizedExperiment::colData(st000284) %>% 
    as.data.frame() %>% 
    tibble::rownames_to_column() 
  testimput <- PomaSummarizedExperiment(features = data, target = target)

  a <- ncol(t(SummarizedExperiment::assay(PomaImpute(testimput, method = "knn", ZerosAsNA = FALSE, RemoveNA = FALSE, cutoff = 8))))
  b <- ncol(t(SummarizedExperiment::assay(PomaImpute(testimput, method = "knn", ZerosAsNA = TRUE, RemoveNA = FALSE, cutoff = 8))))
  c <- ncol(t(SummarizedExperiment::assay(PomaImpute(testimput, method = "knn", ZerosAsNA = FALSE, RemoveNA = TRUE, cutoff = 8))))
  d <- ncol(t(SummarizedExperiment::assay(PomaImpute(testimput, method = "knn", ZerosAsNA = TRUE, RemoveNA = TRUE, cutoff = 8))))

  e <- ncol(t(SummarizedExperiment::assay(PomaImpute(testimput, method = "knn", ZerosAsNA = FALSE, RemoveNA = TRUE, cutoff = 20))))
  f <- ncol(t(SummarizedExperiment::assay(PomaImpute(testimput, method = "knn", ZerosAsNA = FALSE, RemoveNA = TRUE, cutoff = 10))))

  g <- PomaImpute(testimput, method = "half_min", ZerosAsNA = FALSE, RemoveNA = TRUE, cutoff = 20)
  h <- PomaImpute(testimput, method = "knn", ZerosAsNA = FALSE, RemoveNA = TRUE, cutoff = 20)
  
  i <- PomaImpute(testimput, method = "half_min", ZerosAsNA = FALSE, RemoveNA = FALSE, cutoff = 1)
  j <- PomaImpute(testimput, method = "mean", ZerosAsNA = FALSE, RemoveNA = FALSE, cutoff = 1)
  k <- PomaImpute(testimput, method = "median", ZerosAsNA = FALSE, RemoveNA = FALSE, cutoff = 1)

  l <- PomaImpute(testimput, method = "knn", ZerosAsNA = FALSE, RemoveNA = TRUE, cutoff = 20)
  m <- PomaImpute(testimput, method = "knn")

  n <- PomaImpute(testimput, method = "none", RemoveNA = FALSE, cutoff = 2)
  o <- PomaImpute(testimput, method = "none", RemoveNA = FALSE, cutoff = 5)
  p <- PomaImpute(testimput, method = "none", cutoff = 20)
  q <- PomaImpute(testimput, method = "min", cutoff = 20)

  data2 <- t(SummarizedExperiment::assay(testimput))

  data2[1:4, 5] <- 1000
  data2[73:77, 5] <- 1000

  testimput2 <- PomaSummarizedExperiment(features = data2, target = target)

  r <- PomaImpute(testimput2, method = "half_min")
  s <- PomaImpute(testimput2, method = "median")
  t <- PomaImpute(testimput2, method = "mean")
  u <- PomaImpute(testimput2, method = "min")
  v <- PomaImpute(testimput2, method = "knn")
  
  SummarizedExperiment::assay(testimput)[5, 175:190] <- NA
  h_1 <- PomaImpute(testimput, method = "knn", ZerosAsNA = FALSE, RemoveNA = TRUE, cutoff = 1)
  
  ##
  
  expect_equal(testimput@NAMES[1], h@NAMES[1])
  expect_equal(length(testimput@NAMES), length(h@NAMES))
  expect_false(length(testimput@NAMES) == length(h_1@NAMES))
  
  expect_equal(testimput@NAMES[1], n@NAMES[1])
  expect_equal(length(testimput@NAMES), length(n@NAMES))
  
  ##

  expect_equal(a, b)
  expect_equal(b, c)
  expect_equal(a, c)
  expect_equal(a, d)
  expect_equal(b, d)
  expect_equal(c, d)

  expect_equal(d, e)
  expect_equal(c, e)
  expect_equal(e, f)
  
  expect_false(all(assay(g) == assay(h)))
  expect_equal(dim(g), dim(h))

  expect_equal(dim(i), dim(j))
  expect_equal(dim(j), dim(k))

  expect_false(all(assay(i) == assay(j)))
  expect_false(all(assay(j) == assay(k)))
  expect_false(all(assay(k) == assay(i)))

  expect_equal(SummarizedExperiment::assay(l), SummarizedExperiment::assay(m))
  expect_equal(SummarizedExperiment::assay(n), SummarizedExperiment::assay(o))
  expect_true(all(assay(p) == assay(q)))

  expect_equal(dim(r), dim(s))
  expect_equal(dim(s), dim(t))
  expect_equal(dim(t), dim(u))
  expect_equal(dim(u), dim(v))

  ####

  expect_false(all(assay(r) == assay(s)))
  expect_false(all(assay(r) == assay(t)))
  expect_false(all(assay(r) == assay(u)))
  expect_false(all(assay(r) == assay(v)))

  expect_false(all(assay(s) == assay(t)))
  expect_false(all(assay(s) == assay(u)))
  expect_false(all(assay(s) == assay(v)))

  expect_false(all(assay(t) == assay(u)))
  expect_false(all(assay(t) == assay(v)))

  expect_false(all(assay(u) == assay(v)))

  ####

  expect_error(PomaImpute(testimput, method = "non"))
  expect_message(PomaImpute(testimput))
  expect_message(PomaImpute(testimput2))

  ####
  
  expect_message(PomaImpute(st000284, method = "knn"))
  expect_message(PomaImpute(st000284, method = "rf"))
  
  ##
  
  expect_error(PomaImpute(method = "knn"))
  expect_error(PomaImpute(iris, method = "knn"))
  
})

##################################################################
##################################################################

# rfImpute fails many times in virtual machines because it generates a 
# huge proximity matrix that sometimes needs >2 cores to run
# 
# test_that("PomaImpute works skip on Appveyor", {
#   
#   skip_on_appveyor() # rfImpute needs more than 2 cores to run and Appveyor only have 2
#   
#   data("st000336")
#   
#   a_2 <- PomaImpute(st000336, method = "half_min")
#   b_2 <- PomaImpute(st000336, method = "median")
#   c_2 <- PomaImpute(st000336, method = "mean")
#   d_2 <- PomaImpute(st000336, method = "min")
#   e_2 <- PomaImpute(st000336, method = "knn")
#   f_2 <- PomaImpute(st000336, method = "rf")
#   
#   ##
#   
#   expect_false(all(SummarizedExperiment::assay(a_2) == SummarizedExperiment::assay(b_2)))
#   expect_false(all(SummarizedExperiment::assay(b_2) == SummarizedExperiment::assay(c_2)))
#   expect_false(all(SummarizedExperiment::assay(c_2) == SummarizedExperiment::assay(d_2)))
#   expect_false(all(SummarizedExperiment::assay(d_2) == SummarizedExperiment::assay(e_2)))
#   expect_false(all(SummarizedExperiment::assay(e_2) == SummarizedExperiment::assay(f_2)))
#   expect_false(all(SummarizedExperiment::assay(f_2) == SummarizedExperiment::assay(a_2)))
#   
#   ##
#   
#   expect_error(PomaImpute(st000284, method = "rf"))
#   
# })

