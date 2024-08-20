
#' Lasso, Ridge, and Elasticnet Regularized Generalized Linear Models for Binary Outcomes
#'
#' @description `PomaLasso` performs LASSO, Ridge, and Elasticnet regression for feature selection and prediction purposes for binary outcomes.
#'
#' @param data A `SummarizedExperiment` object.
#' @param alpha Numeric. Indicates the elasticnet mixing parameter. alpha = 1 is the LASSO penalty and alpha = 0 is the Ridge penalty.
#' @param ntest Numeric. Indicates the percentage of observations that will be used as test set. Default is NULL (no test set).
#' @param nfolds Numeric. Indicates number of folds for cross-validation (default is 10). Although nfolds can be as large as the sample size (leave-one-out CV), it is not recommended for large datasets. Smallest value allowable is nfolds = 3.
#' @param lambda Numeric. Indicates the user supplied lambda sequence. Typical usage is to have the program compute its own lambda sequence based on `nlambda` and `lambda.min.ratio`. See `?glmnet::glmnet()`.
#' @param labels Logical. Indicates if feature names should be plotted in coefficient plot or not. Default is FALSE.
#' 
#' @export
#'
#' @return A `list` with results.
#' @references Jerome Friedman, Trevor Hastie, Robert Tibshirani (2010). Regularization Paths for Generalized Linear Models via Coordinate Descent. Journal of Statistical Software, 33(1), 1-22. URL http://www.jstatsoft.org/v33/i01/.
#' @author Pol Castellano-Escuder
#'
#' @importFrom magrittr %>%
#' 
#' @examples 
#' data <- POMA::st000336 %>% # Example SummarizedExperiment object included in POMA
#'   PomaImpute() %>% 
#'   PomaNorm()
#' 
#' ## Output is a list with objects `coefficients` (tibble), `coefficients_plot` (ggplot2 object), `cv_plot` (ggplot2 object), and `model` (cv.glmnet object)
#' # LASSO
#' data %>%
#'   PomaLasso(alpha = 1,
#'             ntest = NULL,
#'             nfolds = 10,
#'             lambda = NULL,
#'             labels = TRUE)
#' 
#' # Elasticnet
#' data %>%
#'   PomaLasso(alpha = 0.5,
#'             ntest = NULL,
#'             nfolds = 10,
#'             lambda = NULL,
#'             labels = TRUE)
#' 
#' # Ridge Regression
#' data %>%
#'   PomaLasso(alpha = 0,
#'             ntest = NULL,
#'             nfolds = 10,
#'             lambda = NULL,
#'             labels = FALSE)
PomaLasso <- function(data,
                      alpha = 1,
                      ntest = NULL,
                      nfolds = 10,
                      lambda = NULL,
                      labels = FALSE){

  if (!is(data, "SummarizedExperiment")){
    stop("data is not a SummarizedExperiment object. \nSee POMA::PomaCreateObject or SummarizedExperiment::SummarizedExperiment")
  }
  if (alpha > 1 | alpha < 0) {
    stop("alpha must be a number between 0 and 1")
  }
  if(!is.null(ntest)){
    if (ntest > 50 | ntest < 5) {
      stop("ntest must be a number between 5 and 50 (%)")
    }
  }
  
  features <- t(SummarizedExperiment::assay(data))
  group_factor <- as.factor(SummarizedExperiment::colData(data)[,1])
  to_lasso <- cbind(group_factor, features)
  
  if (length(table(group_factor)[table(group_factor) != 0]) != 2) {
    stop("Grouping factor must have exactly 2 levels (first column of the metadata file)")
  }

  if (!is.null(ntest)){
    
    repeat {

      idx_test <- sample(1:nrow(to_lasso), (ntest/100) * nrow(to_lasso), replace = FALSE)
      
      test <- to_lasso[idx_test ,]
      test_x <- test[,-1]
      test_y <- test[,1]

      train <- to_lasso[-idx_test ,]
      train_x <- train[,-1]
      train_y <- train[,1]
      
      if(length(levels(as.factor(train_y))) == 2 & length(levels(as.factor(test_y))) == 2){
        break
      }
    }
    
    cv_fit <- glmnet::cv.glmnet(data.matrix(train_x), 
                                as.matrix(train_y), 
                                family = "binomial", 
                                nfolds = nfolds, 
                                lambda = lambda, 
                                alpha = alpha)
    
  } else {
    cv_fit <- glmnet::cv.glmnet(features,
                                group_factor, 
                                family = "binomial", 
                                nfolds = nfolds, 
                                lambda = lambda, 
                                alpha = alpha)
  }

  tidied_cv <- broom::tidy(cv_fit)
  glance_cv <- broom::glance(cv_fit)

  cvlasso <- ggplot2::ggplot(tidied_cv, ggplot2::aes(lambda, estimate)) +
    ggplot2::geom_line(color = "blue") +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = conf.low, ymax = conf.high), alpha = 0.2) +
    ggplot2::scale_x_log10() +
    ggplot2::labs(x = "log10(Lambda)",
                  y = "Estimate") +
    ggplot2::geom_vline(xintercept = glance_cv$lambda.min, lty = 2) +
    theme_poma()

  tmp_coeffs <- glmnet::coef.glmnet(cv_fit, s = "lambda.min")
  final_coef <- data.frame(feature = tmp_coeffs@Dimnames[[1]][tmp_coeffs@i + 1], coefficient = tmp_coeffs@x) %>% 
    dplyr::as_tibble()

  if (!is.null(ntest)){
    lasso_pred <- predict(cv_fit, s = cv_fit$lambda.min, newx = data.matrix(test_x), type = "class")
    cm <- caret::confusionMatrix(as.factor(lasso_pred), as.factor(test_y))
  }
  
  tidied_cv2 <- broom::tidy(cv_fit$glmnet.fit)
  tidied_cv2_names <- tidied_cv2 %>% 
    dplyr::arrange(dplyr::desc(abs(estimate))) %>% 
    dplyr::group_by(term) %>% 
    dplyr::slice(1)
  
  coefficientplot <- ggplot2::ggplot(tidied_cv2, ggplot2::aes(lambda, estimate, color = term)) +
    ggplot2::scale_x_log10() +
    ggplot2::labs(x = "log10(Lambda)",
                  y = "Coefficients") +
    ggplot2::geom_line() +
    ggplot2::geom_vline(xintercept = glance_cv$lambda.min, lty = 2) +
    ggplot2::theme_bw() +
    {if(labels)ggplot2::geom_label(data = tidied_cv2_names, ggplot2::aes(label = term))} +
    theme_poma(legend_position = "none") +
    scale_color_poma_d()

  if(!is.null(ntest)){
    return(list(coefficients = final_coef, 
                coefficients_plot = coefficientplot, 
                cv_plot = cvlasso,
                confusion_matrix = cm,
                train_x = train_x,
                train_y = train_y,
                test_x = test_x,
                test_y = test_y,
                model = cv_fit))
  } else {
    return(list(coefficients = final_coef, 
                coefficients_plot = coefficientplot, 
                cv_plot = cvlasso,
                model = cv_fit))
  }
}

