
#' Lasso and Ridge Regularized Generalized Linear Models for Binary Outcomes
#'
#' @description PomaLasso() is an implementation of the lasso an ridge regression from `glmnet` package for binary outcomes.
#'
#' @param data A MSnSet object. First `pData` column must be the subject group/type.
#' @param method Choose between lasso and ridge regression ("lasso", "ridge").
#' @param nfolds Number of folds for CV (default is 10). Although nfolds can be as large as the sample size (leave-one-out CV), it is not recommended for large datasets. Smallest value allowable is nfolds = 3.
#' @param lambda A user supplied lambda sequence. Typical usage is to have the program compute its own lambda sequence based on `nlambda` and `lambda.min.ratio`. See `?glmnet::glmnet()`.
#'
#' @export
#'
#' @return A list with all results including plots and data frames.
#' @references Jerome Friedman, Trevor Hastie, Robert Tibshirani (2010). Regularization Paths for Generalized Linear Models via Coordinate Descent. Journal of Statistical Software, 33(1), 1-22. URL http://www.jstatsoft.org/v33/i01/.
#' @author Pol Castellano-Escuder
#'
#' @import ggplot2
#' @importFrom broom tidy glance
#' @importFrom glmnet cv.glmnet
#' @importFrom crayon red
#' @importFrom clisymbols symbol
#' @importFrom Biobase varLabels pData exprs
PomaLasso <- function(data,
                      method = "lasso",
                      nfolds = 10,
                      lambda = NULL){

  if (missing(data)) {
    stop(crayon::red(clisymbols::symbol$cross, "data argument is empty!"))
  }
  if(!(class(data) == "MSnSet")){
    stop(paste0(crayon::red(clisymbols::symbol$cross, "data is not a MSnSet object."), 
                " \nSee POMA::PomaMSnSetClass or MSnbase::MSnSet"))
  }
  if (missing(method)) {
    warning("method argument is empty! lasso will be used")
  }
  if (!(method %in% c("lasso", "ridge"))) {
    stop(crayon::red(clisymbols::symbol$cross, "Incorrect value for method argument!"))
  }

  Biobase::varLabels(data)[1] <- "Group"

  if (length(levels(as.factor(Biobase::pData(data)$Group))) > 2) {
    stop(crayon::red(clisymbols::symbol$cross, "You data have more than two groups!"))
  }

  X <- t(Biobase::exprs(data))
  Y <- as.factor(Biobase::pData(data)$Group)

  if (method == "lasso"){

    cv_fit <- cv.glmnet(X, Y, family = "binomial", nfolds = nfolds, lambda = lambda)
  }

  if (method == "ridge"){

    cv_fit <- cv.glmnet(X, Y, family = "binomial", nfolds = nfolds, lambda = lambda, alpha = 0)
  }

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

  tmp_coeffs <- coef(cv_fit, s = "lambda.min")
  final_coef <- data.frame(name = tmp_coeffs@Dimnames[[1]][tmp_coeffs@i + 1], coefficient = tmp_coeffs@x)

  ####

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

  return(list(coefficients = final_coef, coefficientPlot = coefficientplot, cvLassoPlot = cvlasso))

}

