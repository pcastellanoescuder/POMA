context("PomaUnivariate")

test_that("PomaUnivariate works", {

  library(SummarizedExperiment)
  
  data("st000284")
  data("st000336")
  
  # st000284_sub <- st000284[,st000284@colData$factors %in% c("CRC", "Healthy")]
  st000336 <- POMA::PomaImpute(st000336, method = "knn")
  
  dims_for_ttest_and_mann <- c(ncol(t(SummarizedExperiment::assay(st000336))), 9)
  dims_for_aov <- c(ncol(t(SummarizedExperiment::assay(st000284))), 
                    length(levels(as.factor(SummarizedExperiment::colData(st000284)[,1]))) + 6)
  dims_for_krusk <- c(ncol(t(SummarizedExperiment::assay(st000284))), 
                      length(levels(as.factor(SummarizedExperiment::colData(st000284)[,1]))) + 7)
  
  univ_ttest <- PomaUnivariate(st000336, method = "ttest", adjust = "fdr")
  univ_aov <- PomaUnivariate(st000284, covariates = FALSE, method = "anova", adjust = "fdr")
  univ_aov_cov <- PomaUnivariate(st000284, covariates = TRUE, method = "anova", adjust = "fdr")
  univ_aov_cov2 <- PomaUnivariate(st000284, covariates = TRUE, method = "anova", adjust = "bonferroni")

  univ_mann <- PomaUnivariate(st000336, method = "mann", adjust = "fdr")
  univ_kruskal <- PomaUnivariate(st000284, method = "kruskal", adjust = "fdr")
  univ_kruskal2 <- PomaUnivariate(st000284, method = "kruskal", adjust = "BY")

  one_cov1 <- PomaUnivariate(st000336, covariates = FALSE, method = "anova", adjust = "fdr")
  one_cov2 <- PomaUnivariate(st000336, covariates = TRUE, method = "anova", adjust = "fdr")
  
  ##

  expect_equal(dims_for_ttest_and_mann, dim(univ_ttest))
  expect_false(dims_for_aov[2] == dim(univ_aov_cov)[2])
  expect_false(dims_for_aov[2] == dim(univ_aov_cov2)[2])
  expect_equal(dims_for_aov[1], dim(univ_aov_cov)[1])
  expect_equal(dims_for_aov[1], dim(univ_aov_cov2)[1])
  expect_equal(dims_for_aov, dim(univ_aov))
  expect_false(all(univ_aov_cov == univ_aov_cov2))

  ##

  expect_equal(dims_for_ttest_and_mann, dim(univ_mann))
  expect_equal(dim(univ_kruskal), dim(univ_kruskal2))
  expect_false(all(univ_kruskal == univ_kruskal2))
  expect_equal(dims_for_krusk, dim(univ_kruskal))
  expect_equal(dims_for_krusk, dim(univ_kruskal2))
  
  expect_false(dim(one_cov1)[2] == dim(one_cov2)[2])
  expect_equal(dim(one_cov1)[1], dim(one_cov2)[1])
  
  ##
  
  expect_error(PomaUnivariate(st000284, covariates = TRUE, method = "anov", adjust = "fdr"))
  expect_error(PomaUnivariate(st000284, covariates = TRUE, adjust = "fdr"))
  
  expect_error(PomaUnivariate(st000336, method = "ttest", adjust = "fd"))
  
  SummarizedExperiment::colData(st000284) <- SummarizedExperiment::colData(st000284)[1]
  expect_error(PomaUnivariate(st000284, method = "anova", covariates = TRUE, adjust = "fdr"))
  
  ##
  
  expect_error(PomaUnivariate(method = "ttest"))
  expect_error(PomaUnivariate(iris, method = "ttest"))
  
})

