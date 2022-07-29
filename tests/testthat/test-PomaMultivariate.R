context("PomaMultivariate")

test_that("PomaMultivariate works", {

  data("st000284")

  #### PCA

  multivariate_pca_1 <- PomaMultivariate(st000284, method = "pca", components = 4,
                                    center = FALSE, scale = FALSE, labels = TRUE)

  multivariate_pca_2 <- PomaMultivariate(st000284, method = "pca", components = 5,
                                    center = TRUE, scale = TRUE)

  ##

  expect_equal(nrow(multivariate_pca_1$scores), nrow(multivariate_pca_2$scores))
  expect_false(ncol(multivariate_pca_1$scores) == ncol(multivariate_pca_2$scores))
  expect_equal(ncol(multivariate_pca_1$scores), 4)
  expect_equal(ncol(multivariate_pca_2$scores), 5)
  
  expect_equal(ncol(multivariate_pca_1$eigenvalues), ncol(multivariate_pca_2$eigenvalues))
  expect_false(nrow(multivariate_pca_1$eigenvalues) == nrow(multivariate_pca_2$eigenvalues))

  ##
  
  df_a <- ggplot2::layer_data(multivariate_pca_1$screeplot)
  df_b <- ggplot2::layer_data(multivariate_pca_1$scoresplot)
  
  df_c <- ggplot2::layer_data(multivariate_pca_2$screeplot)
  df_d <- ggplot2::layer_data(multivariate_pca_2$scoresplot)

  expect_false(length(df_a$y) == length(df_c$y))
  expect_false(length(df_b$y) == length(df_d$y))
  
  #### PLSDA

  multivariate_plsda_1 <- PomaMultivariate(st000284, method = "plsda", components = 3,
                                           center = TRUE, scale = TRUE,
                                           validation = "Mfold", folds = 5, nrepeat = 10, labels = TRUE)

  multivariate_plsda_2 <- PomaMultivariate(st000284, method = "plsda", components = 4,
                                           center = TRUE, scale = TRUE,
                                           validation = "loo", folds = 5, nrepeat = 1, vip = 1)

  ##

  expect_equal(ncol(multivariate_plsda_1$errors_plsda), ncol(multivariate_plsda_2$errors_plsda))
  expect_false(nrow(multivariate_plsda_1$errors_plsda) == nrow(multivariate_plsda_2$errors_plsda))

  expect_false(ncol(multivariate_plsda_1$vip_plsda) == ncol(multivariate_plsda_2$vip_plsda))
  expect_equal(nrow(multivariate_plsda_1$vip_plsda), nrow(multivariate_plsda_2$vip_plsda))
  
  expect_false(ncol(multivariate_plsda_1$scores) == ncol(multivariate_plsda_2$scores))

  ##
  
  df_a <- ggplot2::layer_data(multivariate_plsda_1$scoresplot)
  df_b <- ggplot2::layer_data(multivariate_plsda_1$errors_plsda_plot)
  df_c <- ggplot2::layer_data(multivariate_plsda_1$vip_plsda_plot)

  df_d <- ggplot2::layer_data(multivariate_plsda_2$scoresplot)
  df_e <- ggplot2::layer_data(multivariate_plsda_2$errors_plsda_plot)
  df_f <- ggplot2::layer_data(multivariate_plsda_2$vip_plsda_plot)

  expect_false(ncol(df_a) == ncol(df_d))
  expect_equal(ncol(df_b$y), ncol(df_e$y))
  expect_false(length(df_c$y) == length(df_f$y))

  #### SPLSDA

  multivariate_splsda_1 <- PomaMultivariate(st000284, method = "splsda", components = 3,
                                            center = TRUE, scale = TRUE,
                                            validation = "Mfold", folds = 5, nrepeat = 10,
                                            num_features = 10, labels = TRUE)

  multivariate_splsda_2 <- PomaMultivariate(st000284, method = "splsda", components = 4,
                                            center = TRUE, scale = TRUE,
                                            validation = "Mfold", folds = 5, nrepeat = 10,
                                            num_features = 5)

  ##

  expect_false(nrow(multivariate_splsda_1$selected_variables) ==
                 nrow(multivariate_splsda_2$selected_variables))
  expect_false(nrow(multivariate_splsda_1$errors_splsda) ==
                 nrow(multivariate_splsda_2$errors_splsda))

  expect_true(is.numeric(multivariate_splsda_1$ncomp))
  expect_true(is.numeric(multivariate_splsda_2$ncomp))

  df_a <- ggplot2::layer_data(multivariate_splsda_1$errors_splsda_plot)
  df_b <- ggplot2::layer_data(multivariate_splsda_2$errors_splsda_plot)
  df_c <- ggplot2::layer_data(multivariate_splsda_1$scoresplot)
  df_d <- ggplot2::layer_data(multivariate_splsda_2$scoresplot)

  expect_false(nrow(df_a) == nrow(df_b))
  expect_false(length(df_c$y) == length(df_d$y))
  
  ## ERRORS AND WARNINGS
  
  expect_error(PomaMultivariate(method = "splsda"))
  expect_error(PomaMultivariate(iris, method = "splsda"))
  expect_error(PomaMultivariate(st000284, method = "pc", components = 5))
  expect_error(PomaMultivariate(st000284))
  expect_error(PomaMultivariate(st000284, method = "plsda", validation = "Mfo"))
  expect_message(PomaMultivariate(st000284, method = "plsda"))
  expect_error(PomaMultivariate(st000284, method = "pca", load_length = 2.1))
  expect_error(PomaMultivariate(st000284, method = "pca", load_length = 0.9))
  
})

