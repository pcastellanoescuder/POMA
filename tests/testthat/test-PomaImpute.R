context("PomaImpute")

test_that("PomaImpute works", {

  data <- vroom::vroom("data_ST000284/MET_CRC_ST000284.csv", delim = ",")
  data[,3:ncol(data)] <- data[,3:ncol(data)]*round(runif(n=1, min = 0.01, max = 0.99), 3) # just to create decimals
  data[1:4, 5] <- 0 # create some zeros in one group
  data[73:77, 5] <- 0 # create some zeros in the other group

  data[10:14, 5] <- NA # create some NA in one group (5/66 = 7.6% of NA)
  data[78:81, 5] <- NA # create some NA in the other group (4/66 = 6.1% of NA)

  a <- ncol(PomaImpute(data, method = "knn", ZerosAsNA = F, RemoveNA = F, cutoff = 8))
  b <- ncol(PomaImpute(data, method = "knn", ZerosAsNA = T, RemoveNA = F, cutoff = 8))
  c <- ncol(PomaImpute(data, method = "knn", ZerosAsNA = F, RemoveNA = T, cutoff = 8))
  d <- ncol(PomaImpute(data, method = "knn", ZerosAsNA = T, RemoveNA = T, cutoff = 8))

  e <- ncol(PomaImpute(data, method = "knn", ZerosAsNA = F, RemoveNA = T, cutoff = 20))
  f <- ncol(PomaImpute(data, method = "knn", ZerosAsNA = F, RemoveNA = T, cutoff = 10))

  g <- PomaImpute(data, method = "half_min", ZerosAsNA = F, RemoveNA = T, cutoff = 20)
  h <- PomaImpute(data, method = "knn", ZerosAsNA = F, RemoveNA = T, cutoff = 20)

  i <- PomaImpute(data, method = "half_min", ZerosAsNA = F, RemoveNA = F, cutoff = 1)
  j <- PomaImpute(data, method = "mean", ZerosAsNA = F, RemoveNA = F, cutoff = 1)
  k <- PomaImpute(data, method = "median", ZerosAsNA = F, RemoveNA = F, cutoff = 1)

  l <- PomaImpute(data, method = "knn", ZerosAsNA = F, RemoveNA = T, cutoff = 20)
  m <- PomaImpute(data)

  n <- PomaImpute(data, method = "none", RemoveNA = F, cutoff = 2)
  o <- PomaImpute(data, method = "none", RemoveNA = F, cutoff = 5)
  p <- PomaImpute(data, method = "none", cutoff = 20)
  q <- PomaImpute(data, method = "min", cutoff = 20)

  data[1:4, 5] <- 1000
  data[73:77, 5] <- 1000

  r <- PomaImpute(data, method = "half_min")
  s <- PomaImpute(data, method = "median")
  t <- PomaImpute(data, method = "mean")
  u <- PomaImpute(data, method = "min")
  v <- PomaImpute(data, method = "knn")

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

  expect_false(all(g == h))
  expect_equal(dim(g), dim(h))

  expect_equal(dim(i), dim(j))
  expect_equal(dim(j), dim(k))

  expect_false(all(i == j))
  expect_false(all(j == k))
  expect_false(all(k == i))

  expect_equal(l, m)
  expect_equal(n, o)
  expect_equal(p, q)

  expect_equal(dim(r), dim(s))
  expect_equal(dim(s), dim(t))
  expect_equal(dim(t), dim(u))
  expect_equal(dim(u), dim(v))

  ####

  expect_false(all(r == s))
  expect_false(all(r == t))
  expect_false(all(r == u))
  expect_false(all(r == v))

  expect_false(all(s == t))
  expect_false(all(s == u))
  expect_false(all(s == v))

  expect_false(all(t == u))
  expect_false(all(t == v))

  expect_false(all(u == v))

  ####

  expect_error(PomaImpute(data, method = "non"))
  expect_warning(PomaImpute(data))

})

