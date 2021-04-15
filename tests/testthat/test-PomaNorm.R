context("PomaNorm")

test_that("PomaNorm works", {

  data("st000284")

  data <- t(MSnbase::exprs(st000284))

  data <- data*round(runif(n=1, min = 0.01, max = 0.99), 3) # just to create decimals
  data[1:4, 5] <- 0 # create some zeros in one group
  data[73:77, 5] <- 0 # create some zeros in the other group

  data[10:14, 5] <- NA # create some NA in one group (5/66 = 7.6% of NA)
  data[78:81, 5] <- NA # create some NA in the other group (4/66 = 6.1% of NA)

  data[,1] <- 0 # create column of only zeros
  data[,2] <- 100 # create feature with var = 0

  target <- pData(st000284) %>% rownames_to_column() %>% as.data.frame()
  testnorm <- PomaMSnSetClass(features = data, target = target)

  newdata <- POMA::PomaImpute(testnorm, method = "knn", ZerosAsNA = FALSE, RemoveNA = TRUE, cutoff = 2)
  newdata2 <- POMA::PomaNorm(newdata, method = "log_pareto")

  ####

  a <- dim(PomaNorm(newdata, method = "auto_scaling", round = 2))
  b <- dim(PomaNorm(newdata, method = "level_scaling", round = 2))
  c <- dim(PomaNorm(newdata, method = "log_scaling", round = 2))
  d <- dim(PomaNorm(newdata, method = "log_transformation", round = 2))
  e <- dim(PomaNorm(newdata, method = "vast_scaling", round = 2))
  f <- dim(PomaNorm(newdata, method = "log_pareto", round = 2))
  g <- dim(PomaNorm(newdata, method = "none", round = 2))

  expect_equal(a, b)
  expect_equal(b, c)
  expect_equal(c, d)
  expect_equal(e, f)
  expect_equal(f, g)

  expect_error(PomaNorm(newdata, method = "log", round = 2))
  expect_message(PomaNorm(newdata))

  expect_true(newdata2@processingData@cleaned)
  expect_true(newdata2@processingData@normalised)

  ##
  
  expect_error(PomaNorm(method = "auto_scaling"))
  expect_error(PomaNorm(iris, method = "auto_scaling"))
  
})

