context("PomaImpute")

test_that("PomaImpute works", {

  data <- vroom::vroom("data_ST000284/MET_CRC_ST000284.csv", delim = ",")
  data[,3:ncol(data)] <- data[,3:ncol(data)]*round(runif(n=1, min = 0.01, max = 0.99), 3) # just to create decimals
  data[1:4, 5] <- 0 # create some zeros in one group
  data[73:77, 5] <- 0 # create some zeros in the other group

  a <- ncol(PomaImpute(data, method = "knn", ZerosAsNA = F, RemoveNA = F, cutoff = 2))
  b <- ncol(PomaImpute(data, method = "knn", ZerosAsNA = T, RemoveNA = F, cutoff = 2))
  c <- ncol(PomaImpute(data, method = "knn", ZerosAsNA = F, RemoveNA = T, cutoff = 2))
  d <- ncol(PomaImpute(data, method = "knn", ZerosAsNA = T, RemoveNA = T, cutoff = 2))
  e <- ncol(PomaImpute(data, method = "knn", ZerosAsNA = T, RemoveNA = T, cutoff = 20))
  f <- ncol(PomaImpute(data, method = "knn", ZerosAsNA = T, RemoveNA = T, cutoff = 10))

  g <- PomaImpute(data, method = "half_min", ZerosAsNA = T, RemoveNA = T, cutoff = 20)
  h <- PomaImpute(data, method = "knn", ZerosAsNA = T, RemoveNA = T, cutoff = 20)

  i <- PomaImpute(data, method = "half_min", ZerosAsNA = T, RemoveNA = F, cutoff = 1)
  j <- PomaImpute(data, method = "mean", ZerosAsNA = T, RemoveNA = F, cutoff = 1)
  k <- PomaImpute(data, method = "median", ZerosAsNA = T, RemoveNA = F, cutoff = 1)

  l <- PomaImpute(data, method = "knn", ZerosAsNA = T, RemoveNA = T, cutoff = 20)
  m <- PomaImpute(data)

  n <- PomaImpute(data, method = "none", ZerosAsNA = T, RemoveNA = T, cutoff = 2)
  o <- PomaImpute(data, method = "none", cutoff = 1)
  p <- PomaImpute(data, method = "none", cutoff = 20)
  q <- PomaImpute(data, method = "min", cutoff = 20)

  ####

  expect_equal(a, b)
  expect_equal(b, c)

  expect_true(c != d)
  expect_true(d != e)

  expect_equal(c, e)
  expect_equal(e, f)

  expect_false(all(g == h))
  expect_equal(dim(g), dim(h))

  expect_equal(dim(i), dim(j))
  expect_equal(dim(j), dim(k))

  expect_false(all(i == j))
  expect_false(all(j == k))
  expect_false(all(k == i))

  expect_equal(l, m)
  expect_equal(n, o)

  expect_false(all(p == q))
  expect_equal(dim(p), dim(q))

  expect_error(PomaImpute(data, method = "non"))
  expect_warning(PomaImpute(data))

})

