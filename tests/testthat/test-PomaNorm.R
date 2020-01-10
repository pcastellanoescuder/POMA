context("PomaNorm")

test_that("PomaNorm works", {

  data("st000284")

  a <- dim(PomaNorm(st000284, method = "auto_scaling", round = 2))
  b <- dim(PomaNorm(st000284, method = "level_scaling", round = 2))
  c <- dim(PomaNorm(st000284, method = "log_scaling", round = 2))
  d <- dim(PomaNorm(st000284, method = "log_transformation", round = 2))
  e <- dim(PomaNorm(st000284, method = "vast_scaling", round = 2))
  f <- dim(PomaNorm(st000284, method = "log_pareto", round = 2))
  g <- dim(PomaNorm(st000284, method = "none", round = 2))

  expect_equal(a, b)
  expect_equal(b, c)
  expect_equal(c, d)
  expect_equal(e, f)
  expect_equal(f, g)

  expect_error(PomaNorm(st000284, method = "log", round = 2))
  expect_error(PomaNorm(st000284))

})

