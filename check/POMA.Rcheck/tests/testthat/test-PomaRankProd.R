context("PomaRankProd")

test_that("PomaRankProd works", {

  data("st000284")
  
  target <- Biobase::pData(st000284)[1:100,] %>% tibble::rownames_to_column("ID")
  e <- Biobase::exprs(st000284)[,1:100]
  
  data <- PomaMSnSetClass(target = target, features = t(e))
  
  toy_data <- POMA::PomaNorm(data, method = "log_scaling")
  Biobase::pData(toy_data)$groups <- c(rep("C", 25), rep("G", 25))

  ##

  RP_one <- PomaRankProd(data, logged = TRUE, logbase = 2)
  RP_two <- PomaRankProd(data, logged = TRUE, logbase = 10)

  RP_five <- PomaRankProd(data, cutoff = 0.05, method = "pfp")
  RP_six <- PomaRankProd(data, cutoff = 0.05, method = "pval")

  ##

  expect_error(PomaRankProd())
  expect_error(PomaRankProd(data, method = "pfd"))
  expect_error(PomaRankProd(toy_data))
  expect_warning(PomaRankProd(data))
  expect_error(PomaRankProd(iris))

  ##

  expect_equal(dim(RP_one$upregulated), dim(RP_two$upregulated))
  expect_equal(dim(RP_one$downregulated), dim(RP_two$downregulated))

  expect_false(dim(RP_five$upregulated)[1] == dim(RP_six$upregulated)[1])
  expect_false(dim(RP_five$downregulated)[1] == dim(RP_six$downregulated)[1])
  
})

