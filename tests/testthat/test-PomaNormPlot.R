context("PomaNormPlot")

test_that("PomaNormPlot works", {

  library(tidyverse)

  data <- vroom::vroom("data_ST000284/MET_CRC_ST000284.csv", delim = ",")

  norm_none <- PomaNorm(data, method = "none")
  norm_ls <- PomaNorm(data, method = "log_scaling")

  a <- PomaNormPlot(norm_none)
  b <- PomaNormPlot(norm_ls)
  c <- PomaNormPlot(norm_none, group = "metabolites")
  d <- PomaNormPlot(norm_ls, group = "metabolites")

  e <- PomaNormPlot(norm_none, group = "subjects")

  df_a <- layer_data(a)
  df_b <- layer_data(b)
  df_c <- layer_data(c)
  df_d <- layer_data(d)
  df_e <- layer_data(e)

  ####

  expect_false(all(df_a$ymin == df_c$ymin))
  expect_false(all(df_b$ymin == df_d$ymin))
  expect_false(all(df_a$ymin == df_b$ymin))
  expect_false(all(df_c$ymin == df_d$ymin))

  expect_equal(df_a, df_e)

  expect_warning(PomaNormPlot(norm_none))
  expect_warning(PomaNormPlot(norm_ls))

})

