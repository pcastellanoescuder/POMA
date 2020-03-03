context("PomaRankProd")

test_that("PomaRankProd works", {

  data("st000284")

  toy_data <- POMA::PomaNorm(st000284, method = "log_scaling")
  Biobase::pData(toy_data)$groups <- c(rep("C", 100), rep("H", 29), rep("G", 29))

  ##

  RP_one <- PomaRankProd(st000284, logged = TRUE, logbase = 2)
  RP_two <- PomaRankProd(st000284, logged = TRUE, logbase = 10)

  RP_five <- PomaRankProd(st000284, cutoff = 0.05, method = "pfp")
  RP_six <- PomaRankProd(st000284, cutoff = 0.05, method = "pval")

  ##

  expect_error(PomaRankProd())
  expect_error(PomaRankProd(st000284, method = "pfd"))
  expect_error(PomaRankProd(toy_data))
  expect_warning(PomaRankProd(st000284))

  ##

  expect_equal(dim(RP_one$upregulated), dim(RP_two$upregulated))
  expect_equal(dim(RP_one$downregulated), dim(RP_two$downregulated))

  expect_false(dim(RP_five$upregulated)[1] == dim(RP_six$upregulated)[1])
  expect_false(dim(RP_five$downregulated)[1] == dim(RP_six$downregulated)[1])

})

