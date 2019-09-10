context("PomaNorm")

test_that("PomaNorm works", {

  data <- vroom::vroom("data_ST000284/MET_CRC_ST000284.csv", delim = ",")

  a <- dim(PomaNorm(data, method = "auto_scaling", round = 2))
  b <- dim(PomaNorm(data, method = "level_scaling", round = 2))
  c <- dim(PomaNorm(data, method = "log_scaling", round = 2))
  d <- dim(PomaNorm(data, method = "log_transformation", round = 2))
  e <- dim(PomaNorm(data, method = "vast_scaling", round = 2))
  f <- dim(PomaNorm(data, method = "log_pareto", round = 2))
  g <- dim(PomaNorm(data, method = "none", round = 2))

  expect_equal(a, b)
  expect_equal(b, c)
  expect_equal(c, d)
  expect_equal(e, f)
  expect_equal(f, g)

  expect_error(PomaNorm(data, method = "log", round = 2))
  expect_error(PomaNorm(data))

})

