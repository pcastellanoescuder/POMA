context("PomaImpute")

test_that("PomaImpute works", {

  data("st000284")

  data <- t(Biobase::exprs(st000284))

  data <- data*round(runif(n = 1, min = 0.01, max = 0.99), 3) # just to create decimals
  data[1:4, 5] <- 0 # create some zeros in one group
  data[73:77, 5] <- 0 # create some zeros in the other group

  data[10:14, 5] <- NA # create some NA in one group (5/66 = 7.6% of NA)
  data[78:81, 5] <- NA # create some NA in the other group (4/66 = 6.1% of NA)

  target <- pData(st000284) %>% rownames_to_column() %>% as.data.frame()
  testimput <- PomaMSnSetClass(features = data, target = target)

  a <- ncol(t(Biobase::exprs(PomaImpute(testimput, method = "knn", ZerosAsNA = F, RemoveNA = F, cutoff = 8))))
  b <- ncol(t(Biobase::exprs(PomaImpute(testimput, method = "knn", ZerosAsNA = T, RemoveNA = F, cutoff = 8))))
  c <- ncol(t(Biobase::exprs(PomaImpute(testimput, method = "knn", ZerosAsNA = F, RemoveNA = T, cutoff = 8))))
  d <- ncol(t(Biobase::exprs(PomaImpute(testimput, method = "knn", ZerosAsNA = T, RemoveNA = T, cutoff = 8))))

  e <- ncol(t(Biobase::exprs(PomaImpute(testimput, method = "knn", ZerosAsNA = F, RemoveNA = T, cutoff = 20))))
  f <- ncol(t(Biobase::exprs(PomaImpute(testimput, method = "knn", ZerosAsNA = F, RemoveNA = T, cutoff = 10))))

  e_rf <- ncol(t(Biobase::exprs(PomaImpute(testimput, method = "rf", ZerosAsNA = F, RemoveNA = T, cutoff = 20))))
  f_rf <- ncol(t(Biobase::exprs(PomaImpute(testimput, method = "rf", ZerosAsNA = F, RemoveNA = T, cutoff = 10))))
  
  g <- PomaImpute(testimput, method = "half_min", ZerosAsNA = F, RemoveNA = T, cutoff = 20)
  h <- PomaImpute(testimput, method = "knn", ZerosAsNA = F, RemoveNA = T, cutoff = 20)

  i <- PomaImpute(testimput, method = "half_min", ZerosAsNA = F, RemoveNA = F, cutoff = 1)
  j <- PomaImpute(testimput, method = "mean", ZerosAsNA = F, RemoveNA = F, cutoff = 1)
  k <- PomaImpute(testimput, method = "median", ZerosAsNA = F, RemoveNA = F, cutoff = 1)

  l <- PomaImpute(testimput, method = "knn", ZerosAsNA = F, RemoveNA = T, cutoff = 20)
  m <- PomaImpute(testimput, method = "knn")

  n <- PomaImpute(testimput, method = "none", RemoveNA = F, cutoff = 2)
  o <- PomaImpute(testimput, method = "none", RemoveNA = F, cutoff = 5)
  p <- PomaImpute(testimput, method = "none", cutoff = 20)
  q <- PomaImpute(testimput, method = "min", cutoff = 20)

  data2 <- t(Biobase::exprs(testimput))

  data2[1:4, 5] <- 1000
  data2[73:77, 5] <- 1000

  testimput2 <- PomaMSnSetClass(features = data2, target = target)

  r <- PomaImpute(testimput2, method = "half_min")
  s <- PomaImpute(testimput2, method = "median")
  t <- PomaImpute(testimput2, method = "mean")
  u <- PomaImpute(testimput2, method = "min")
  v <- PomaImpute(testimput2, method = "knn")
  y <- PomaImpute(testimput2, method = "rf")
  
  ####

  expect_equal(a, b)
  expect_equal(b, c)
  expect_equal(a, c)
  expect_false(a == d)
  expect_false(b == d)
  expect_false(c == d)

  expect_true(d != e)
  expect_equal(c, e)
  expect_equal(e, f)
  expect_equal(e, e_rf)
  expect_equal(f, f_rf)
  
  expect_false(all(exprs(g) == exprs(h)))
  expect_equal(dim(g), dim(h))

  expect_equal(dim(i), dim(j))
  expect_equal(dim(j), dim(k))

  expect_false(all(exprs(i) == exprs(j)))
  expect_false(all(exprs(j) == exprs(k)))
  expect_false(all(exprs(k) == exprs(i)))

  expect_equal(Biobase::exprs(l), Biobase::exprs(m))
  expect_equal(Biobase::exprs(n), Biobase::exprs(o))
  expect_true(all(exprs(p) == exprs(q)))

  expect_equal(dim(r), dim(s))
  expect_equal(dim(s), dim(t))
  expect_equal(dim(t), dim(u))
  expect_equal(dim(u), dim(v))

  ####

  expect_false(all(exprs(r) == exprs(s)))
  expect_false(all(exprs(r) == exprs(t)))
  expect_false(all(exprs(r) == exprs(u)))
  expect_false(all(exprs(r) == exprs(v)))

  expect_false(all(exprs(s) == exprs(t)))
  expect_false(all(exprs(s) == exprs(u)))
  expect_false(all(exprs(s) == exprs(v)))

  expect_false(all(exprs(t) == exprs(u)))
  expect_false(all(exprs(t) == exprs(v)))

  expect_false(all(exprs(u) == exprs(v)))
  expect_false(all(exprs(v) == exprs(y)))

  ####

  expect_error(PomaImpute(testimput, method = "non"))
  expect_warning(PomaImpute(testimput))
  expect_warning(PomaImpute(testimput2))

  ####

  expect_true(g@processingData@cleaned)
  expect_true(h@processingData@cleaned)

  ####
  
  expect_warning(PomaImpute(st000284, method = "knn"))
  expect_error(PomaImpute(st000284, method = "rf")) # rfImpute doesn't work without NAs
})

