context("PomaVolcano")

test_that("PomaVolcano works", {

  data <- vroom::vroom("data_ST000284/MET_CRC_ST000284.csv", delim = ",")

  univ_ttest <- PomaUnivariate(data, method = "ttest", adjust = "fdr")
  univ_mann <- PomaUnivariate(data, method = "mann", adjust = "fdr")

  a <- PomaVolcano(univ_ttest)
  b <- PomaVolcano(univ_mann, pval = "raw", Pval_cutoff = 0.05, FC_cutoff = 1.5, xlim = 2)
  c <- PomaVolcano(univ_ttest, pval = "raw", Pval_cutoff = 0.05, FC_cutoff = 1.5, xlim = 2)

  df_a <- layer_data(a)
  df_b <- layer_data(b)
  df_c <- layer_data(c)

  expect_false(all(df_a$y == df_b$y))
  expect_equal(df_a, df_c)

  expect_warning(PomaVolcano(univ_ttest))
  expect_error(PomaVolcano(univ_ttest, pval = "ra"))

})
