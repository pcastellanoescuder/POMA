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

  expect_equal(a, b)
  expect_equal(b, c)

  expect_true(c != d)
  expect_true(d != e)

  expect_equal(c, e)
  expect_equal(e, f)

})

