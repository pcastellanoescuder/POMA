
#' Lasso, Ridge and Elasticnet Regularized Generalized Linear Models for Binary Outcomes
#'
#' @description PomaLasso() is an implementation of the lasso, ridge and elasticnet regression from `glmnet` package for binary outcomes.
#'
#' @param data A SummarizedExperiment object.
#' @param alpha Elasticnet mixing parameter. alpha = 1 is the lasso penalty and alpha = 0 is the ridge penalty. This value must be between 0 and 1.
#' @param ntest Numeric indicating the percentage of observations that will be used as test set. Default is NULL (no test set).
#' @param nfolds Number of folds for CV (default is 10). Although nfolds can be as large as the sample size (leave-one-out CV), it is not recommended for large datasets. Smallest value allowable is nfolds = 3.
#' @param lambda A user supplied lambda sequence. Typical usage is to have the program compute its own lambda sequence based on `nlambda` and `lambda.min.ratio`. See `?glmnet::glmnet()`.
#' @param labels Logical indicating if feature names should be plotted in coefficient plot or not. Default is FALSE.
#' 
#' @export
#'
#' @return A list with all results including plots, tables and the resulting prediction model.
#' @references Jerome Friedman, Trevor Hastie, Robert Tibshirani (2010). Regularization Paths for Generalized Linear Models via Coordinate Descent. Journal of Statistical Software, 33(1), 1-22. URL http://www.jstatsoft.org/v33/i01/.
#' @author Pol Castellano-Escuder
#'
#' @importFrom magrittr %>%
#' 
#' @examples 
#' data("st000336")
#' 
#' # lasso
#' st000336 %>%
#'   PomaImpute() %>%
#'   PomaNorm() %>%
#'   PomaOutliers() %>%
#'   PomaLasso()
#' 
#' # elasticnet
#' st000336 %>%
#'   PomaImpute() %>%
#'   PomaNorm() %>%
#'   PomaOutliers() %>%
#'   PomaLasso(alpha = 0.5)
#' 
#' # ridge
#' st000336 %>%
#'   PomaImpute() %>%
#'   PomaNorm() %>%
#'   PomaOutliers() %>%
#'   PomaLasso(alpha = 0)
PomaLasso <- function(data,
                      alpha = 1,
                      ntest = NULL,
                      nfolds = 10,
                      lambda = NULL,
                      labels = FALSE){

  if (missing(data)) {
    stop("data argument is empty!")
  }
  if(!is(data[1], "SummarizedExperiment")){
    stop("data is not a SummarizedExperiment object. \nSee POMA::PomaSummarizedExperiment or SummarizedExperiment::SummarizedExperiment")
  }
  if (alpha > 1 | alpha < 0) {
    stop("alpha must be a number between 0 and 1...")
  }
  if(!is.null(ntest)){
    if (ntest > 50 | ntest < 5) {
      stop("ntest must be a number between 5 and 50...")
    }
  }
  if (length(levels(as.factor(SummarizedExperiment::colData(data)[,1]))) > 2) {
    stop("Your data have more than two groups!")
  }
  if (length(levels(as.factor(SummarizedExperiment::colData(data)[,1]))) < 2) {
    stop("Your data have less than two groups!")
  }

  features <- t(SummarizedExperiment::assay(data))
  response <- as.factor(SummarizedExperiment::colData(data)[,1])
  lasso_data <- cbind(response, features)

  n <- nrow(lasso_data)

  if(!is.null(ntest)){
    
    repeat {

      idx_test <- sample(1:n, (ntest/100)*n, replace = FALSE)
      
      test <- lasso_data[idx_test ,]
      test_x <- test[,-1]
      test_y <- test[,1]

      train <- lasso_data[-idx_test ,]
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
    
  } 
  else {
    cv_fit <- glmnet::cv.glmnet(features,
                                response, 
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
    ggplot2::theme_bw()

  tmp_coeffs <- glmnet::coef.glmnet(cv_fit, s = "lambda.min")
  final_coef <- data.frame(feature = tmp_coeffs@Dimnames[[1]][tmp_coeffs@i + 1], coefficient = tmp_coeffs@x) %>% 
    dplyr::as_tibble()

  if(!is.null(ntest)){
    lasso_pred <- glmnet::predict(cv_fit, s = cv_fit$lambda.min, newx = data.matrix(test_x), type = "class")
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
    ggplot2::theme(legend.position = "none") +
    ggplot2::scale_color_viridis_d(begin = 0, end = 0.8)

  if(!is.null(ntest)){
    return(list(coefficients = final_coef, 
                coefficientPlot = coefficientplot, 
                cvLassoPlot = cvlasso,
                confusionMatrix = cm,
                train_x = train_x,
                train_y = train_y,
                test_x = test_x,
                test_y = test_y,
                model = cv_fit))
  } else {
    return(list(coefficients = final_coef, 
                coefficientPlot = coefficientplot, 
                cvLassoPlot = cvlasso,
                model = cv_fit))
  }

}

