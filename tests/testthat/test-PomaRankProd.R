context("PomaRankProd")

test_that("PomaRankProd works", {

  data("st000284")
  
  target <- SummarizedExperiment::colData(st000284)[1:100,] %>% 
    as.data.frame() %>% 
    tibble::rownames_to_column("ID")
  e <- SummarizedExperiment::assay(st000284)[,1:100]
  
  data <- PomaSummarizedExperiment(target = target, features = t(e))
  
  toy_data <- POMA::PomaNorm(data, method = "log_scaling")
  SummarizedExperiment::colData(toy_data)$groups <- c(rep("C", 25), rep("G", 25))

  ##

  RP_one <- PomaRankProd(data, logged = TRUE, logbase = 2)
  RP_two <- PomaRankProd(data, logged = TRUE, logbase = 10)

  RP_five <- PomaRankProd(data, cutoff = 0.05, method = "pfp")
  RP_six <- PomaRankProd(data, cutoff = 0.05, method = "pval")

  ##

  expect_error(PomaRankProd())
  expect_error(PomaRankProd(data, method = "pfd"))
  expect_error(PomaRankProd(toy_data))
  expect_message(PomaRankProd(data))
  expect_error(PomaRankProd(iris))

  ##

  expect_equal(dim(RP_one$upregulated), dim(RP_two$upregulated))
  expect_equal(dim(RP_one$downregulated), dim(RP_two$downregulated))

  expect_false(dim(RP_five$upregulated)[1] == dim(RP_six$upregulated)[1])
  expect_false(dim(RP_five$downregulated)[1] == dim(RP_six$downregulated)[1])
  
})

