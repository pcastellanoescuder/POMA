context("PomaVolcano")

test_that("PomaVolcano works", {

  data("st000284")

  iris_example <- PomaMSnSetClass(target = data.frame(ID = 1:150, Group = iris$Species), features = iris[,1:4])
    
  a <- PomaVolcano(st000284, pval = "adjusted", adjust = "fdr")
  b <- PomaVolcano(st000284, pval = "adjusted", pval_cutoff = 0.05, log2FC = 0.6, xlim = 2, adjust = "fdr")
  c <- PomaVolcano(st000284, pval = "raw", pval_cutoff = 0.05, log2FC = 0.6, xlim = 2, adjust = "fdr")
  d <- PomaVolcano(st000284, pval = "raw", pval_cutoff = 0.05, log2FC = 0.6, xlim = 2, adjust = "bonferroni")
  
  df_a <- layer_data(a)
  df_b <- layer_data(b)
  df_c <- layer_data(c)
  df_d <- layer_data(d)
  
  ##
  
  expect_equal(df_a, df_b)
  expect_false(all(df_a$y == df_c$y))
  
  expect_equal(df_c$label, df_d$label)

  ##
  
  expect_warning(PomaVolcano(st000284, adjust = "fdr"))
  expect_warning(PomaVolcano(st000284, pval = "raw"))
  
  expect_error(PomaVolcano(st000284, pval = "ra", adjust = "fdr"))
  expect_error(PomaVolcano(st000284, pval = "raw", adjust = "fd"))
  
  expect_error(PomaVolcano(iris_example, pval = "raw", adjust = "fdr"))
  
  ##
  
  expect_error(PomaVolcano())
  expect_error(PomaVolcano(iris))
  
})

