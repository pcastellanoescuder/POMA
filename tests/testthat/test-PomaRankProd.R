context("PomaRankProd")

test_that("PomaRankProd works", {

  library(tidyverse)
  library(RankProd)

  data <- vroom::vroom("data_ST000284/MET_CRC_ST000284.csv", delim = ",")
  data <- PomaImpute(data, method = "knn")

  toy_data <- data.frame(ID = c("One", "Two", "Three"),
                         Group = c("A", "B", "C"),
                         Alanine = c(123, 453, 432))

  ##

  RP_one <- PomaRankProd(data, logged = TRUE, logbase = 2)
  RP_two <- PomaRankProd(data, logged = TRUE, logbase = 10)
  # RP_three <- PomaRankProd(data, logged = FALSE, logbase = 2)
  # RP_four <- PomaRankProd(data, logged = FALSE, logbase = 10)
  #
  # RP_five <- PomaRankProd(data, cutoff = 0.05, method = "pfp")
  # RP_six <- PomaRankProd(data, cutoff = 0.05, method = "pval")
  # RP_seven <- PomaRankProd(data, cutoff = 0.01, method = "pfp")
  # RP_eight <- PomaRankProd(data, cutoff = 0.01, method = "pval")

  ##

  expect_error(PomaRankProd())
  expect_error(PomaRankProd(data, method = "pfd"))
  expect_error(PomaRankProd(toy_data))
  expect_warning(PomaRankProd(data))

  ##

  expect_equal(dim(RP_one$upregulated), dim(RP_two$upregulated))
  expect_equal(dim(RP_one$downregulated), dim(RP_two$downregulated))

})

