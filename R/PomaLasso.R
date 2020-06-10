
#' Lasso, Ridge and Elasticnet Regularized Generalized Linear Models for Binary Outcomes
#'
#' @description PomaLasso() is an implementation of the lasso, ridge and elasticnet regression from `glmnet` package for binary outcomes.
#'
#' @param data A MSnSet object. First `pData` column must be the subject group/type.
#' @param alpha Elasticnet mixing parameter. alpha = 1 is the lasso penalty and alpha = 0 is the ridge penalty. This value must be between 0 and 1.
#' @param ntest Numeric indicating the percentage of observations that will be used as test set. Default is a 20 percent of observations.
#' @param nfolds Number of folds for CV (default is 10). Although nfolds can be as large as the sample size (leave-one-out CV), it is not recommended for large datasets. Smallest value allowable is nfolds = 3.
#' @param lambda A user supplied lambda sequence. Typical usage is to have the program compute its own lambda sequence based on `nlambda` and `lambda.min.ratio`. See `?glmnet::glmnet()`.
#'
#' @export
#'
#' @return A list with all results including plots, data frames and the resulting prediction model.
#' @references Jerome Friedman, Trevor Hastie, Robert Tibshirani (2010). Regularization Paths for Generalized Linear Models via Coordinate Descent. Journal of Statistical Software, 33(1), 1-22. URL http://www.jstatsoft.org/v33/i01/.
#' @author Pol Castellano-Escuder
#'
#' @import ggplot2
#' @importFrom broom tidy glance
#' @importFrom glmnet cv.glmnet
#' @importFrom crayon red
#' @importFrom clisymbols symbol
#' @importFrom caret confusionMatrix
#' @importFrom Biobase varLabels pData exprs
PomaLasso <- function(data,
                      alpha = 1,
                      ntest = 20,
                      nfolds = 10,
                      lambda = NULL){

  if (missing(data)) {
    stop(crayon::red(clisymbols::symbol$cross, "data argument is empty!"))
  }
  if(!(class(data) == "MSnSet")){
    stop(paste0(crayon::red(clisymbols::symbol$cross, "data is not a MSnSet object."), 
                " \nSee POMA::PomaMSnSetClass or MSnbase::MSnSet"))
  }
  if (alpha > 1 | alpha < 0) {
    stop(crayon::red(clisymbols::symbol$cross, "alpha must be a number between 0 and 1..."))
  }
  if (ntest > 50 | ntest < 10) {
    stop(crayon::red(clisymbols::symbol$cross, "ntest must be a number between 10 and 50..."))
  }
  
  Biobase::varLabels(data)[1] <- "Group"

  if (length(levels(as.factor(Biobase::pData(data)$Group))) > 2) {
    stop(crayon::red(clisymbols::symbol$cross, "Your data have more than two groups!"))
  }
  if (length(levels(as.factor(Biobase::pData(data)$Group))) < 2) {
    stop(crayon::red(clisymbols::symbol$cross, "Your data have less than two groups!"))
  }

  features <- t(Biobase::exprs(data))
  response <- as.factor(Biobase::pData(data)$Group)
  lasso_data <- cbind(response, features)

  n <- nrow(lasso_data)
  
  ## MODEL
  if(ntest != 0){
    
    repeat{
      
      ## TEST
      idx_test <- sample(1:n, (ntest/100)*n, replace = FALSE)
      
      test <- lasso_data[idx_test ,]
      test_x <- test[,-1]
      test_y <- test[,1]
      
      ## TRAIN
      train <- lasso_data[-idx_test ,]
      train_x <- train[,-1]
      train_y <- train[,1]
      
      if(length(levels(as.factor(train_y))) == 2 & length(levels(as.factor(test_y)))){
        break
      }
    }
    
    cv_fit <- cv.glmnet(data.matrix(train_x), as.matrix(train_y), family = "binomial", nfolds = nfolds, lambda = lambda, alpha = alpha)
    
  } 
  else {
    cv_fit <- cv.glmnet(features, response, family = "binomial", nfolds = nfolds, lambda = lambda, alpha = alpha)
  }

  ## CROSS-VALIDATION PLOT
  tidied_cv <- broom::tidy(cv_fit)
  glance_cv <- broom::glance(cv_fit)

  cvlasso <- ggplot(tidied_cv, aes(lambda, estimate)) +
    geom_line(color = "blue") +
    geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2) +
    scale_x_log10() +
    xlab("log10(Lambda)") +
    ylab("Estimate") +
    geom_vline(xintercept = glance_cv$lambda.min) +
    geom_vline(xintercept = glance_cv$lambda.1se, lty = 2) +
    theme_bw()

  ## COEFFICIENTS
  tmp_coeffs <- coef(cv_fit, s = "lambda.min")
  final_coef <- data.frame(feature = tmp_coeffs@Dimnames[[1]][tmp_coeffs@i + 1], coefficient = tmp_coeffs@x)

  ## MODEL VALIDATION
  if(ntest != 0){
    lasso_pred <- predict(cv_fit, s = cv_fit$lambda.min, newx = data.matrix(test_x), type = "class")
    cm <- caret::confusionMatrix(as.factor(lasso_pred), as.factor(test_y))
  }
  
  ## COEFFICIENT PLOT
  tidied_cv2 <- broom::tidy(cv_fit$glmnet.fit)

  coefficientplot <- ggplot(tidied_cv2, aes(lambda, estimate, color = term)) +
    scale_x_log10() +
    xlab("log10(Lambda)") +
    ylab("Coefficients") +
    geom_line() +
    geom_vline(xintercept = glance_cv$lambda.min) +
    geom_vline(xintercept = glance_cv$lambda.1se, lty = 2) +
    theme_bw() +
    theme(legend.position = "none")

  if(ntest != 0){
    return(list(coefficients = final_coef, coefficientPlot = coefficientplot, cvLassoPlot = cvlasso,
                confusionMatrix = cm, model = cv_fit))
  } else {
    return(list(coefficients = final_coef, coefficientPlot = coefficientplot, cvLassoPlot = cvlasso,
                model = cv_fit))
  }

}

