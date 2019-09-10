context("PomaMultivariate")

test_that("PomaMultivariate works", {

  library(tidyverse)

  data <- vroom::vroom("data_ST000284/MET_CRC_ST000284.csv", delim = ",")

  #### pca

  multivariate_pca_1 <- PomaMultivariate(data, method = "pca", components = 5,
                                    center = FALSE, scale = FALSE)

  multivariate_pca_2 <- PomaMultivariate(data, method = "pca", components = 5,
                                    center = TRUE, scale = TRUE)

  ####

  expect_equal(dim(multivariate_pca_1$score_data), dim(multivariate_pca_2$score_data))
  expect_false(all(multivariate_pca_1$score_data == multivariate_pca_2$score_data))

  expect_equal(dim(multivariate_pca_1$eigenvalues), dim(multivariate_pca_2$eigenvalues))
  expect_false(all(multivariate_pca_1$eigenvalues == multivariate_pca_2$eigenvalues))

  expect_error(PomaMultivariate(data, method = "pc", components = 5))
  expect_error(PomaMultivariate(data))
  expect_error(PomaMultivariate(data, method = "pca", validation = "Mfo"))
  expect_warning(PomaMultivariate(data, method = "pca"))

  df_a <- layer_data(multivariate_pca_1$screeplot)
  df_b <- layer_data(multivariate_pca_1$scoresplot)
  df_c <- layer_data(multivariate_pca_2$screeplot)
  df_d <- layer_data(multivariate_pca_2$scoresplot)

  expect_false(all(df_a$y == df_c$y))
  expect_false(all(df_b$y == df_d$y))

  #### plsda

  multivariate_plsda_1 <- PomaMultivariate(data, method = "plsda", components = 5,
                                           center = TRUE, scale = TRUE,
                                           validation = "Mfold", folds = 5, nrepeat = 10)

  multivariate_plsda_2 <- PomaMultivariate(data, method = "plsda", components = 4,
                                           center = TRUE, scale = TRUE,
                                           validation = "loo", folds = 5, nrepeat = 1)

  ####

  expect_equal(ncol(multivariate_plsda_1$errors_plsda), ncol(multivariate_plsda_2$errors_plsda))
  expect_false(all(multivariate_plsda_1$errors_plsda$Component ==
                    multivariate_plsda_2$errors_plsda$Component))

  expect_false(ncol(multivariate_plsda_1$plsda_vip_table) == ncol(multivariate_plsda_2$plsda_vip_table))

  expect_false(ncol(multivariate_plsda_1$scores_plsda) == ncol(multivariate_plsda_2$scores_plsda))

  df_a <- layer_data(multivariate_plsda_1$scoresplot)
  df_b <- layer_data(multivariate_plsda_1$errors_plsda_plot)
  df_c <- layer_data(multivariate_plsda_1$vip_plsda_plot)

  df_d <- layer_data(multivariate_plsda_2$scoresplot)
  df_e <- layer_data(multivariate_plsda_2$errors_plsda_plot)
  df_f <- layer_data(multivariate_plsda_2$vip_plsda_plot)

  expect_equal(df_a$y, df_d$y)
  expect_false(all(df_b$y == df_e$y))
  expect_equal(df_c$y, df_f$y)

  #### splsda

  multivariate_splsda_1 <- PomaMultivariate(data, method = "splsda", components = 5,
                                            center = TRUE, scale = TRUE,
                                            validation = "Mfold", folds = 5, nrepeat = 10,
                                            num_features = 10)

  multivariate_splsda_2 <- PomaMultivariate(data, method = "splsda", components = 5,
                                            center = TRUE, scale = TRUE,
                                            validation = "Mfold", folds = 5, nrepeat = 10,
                                            num_features = 5)

  ####

  expect_false(nrow(multivariate_splsda_1$selected_variables) ==
                 nrow(multivariate_splsda_2$selected_variables))
  expect_false(nrow(multivariate_splsda_1$errors_splsda) ==
                 nrow(multivariate_splsda_2$errors_splsda))

  expect_true(is.numeric(multivariate_splsda_1$ncomp))
  expect_true(is.numeric(multivariate_splsda_2$ncomp))

  df_a <- layer_data(multivariate_splsda_1$bal_error_rate)
  df_b <- layer_data(multivariate_splsda_2$bal_error_rate)
  df_c <- layer_data(multivariate_splsda_1$splsda_scores_plot)
  df_d <- layer_data(multivariate_splsda_2$splsda_scores_plot)

  expect_false(all(df_a$y == df_b$y))
  expect_false(all(df_c$y == df_d$y))

})

