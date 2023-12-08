
#' Classification Random Forest
#'
#' @description `PomaRandForest` performs classification random forest. This method can be used both for prediction and variable selection.
#'
#' @param data A `SummarizedExperiment` object.
#' @param ntest Numeric. Indicates the percentage of observations that will be used as test set. Default is NULL (no test set).
#' @param ntree Numeric. Indicates the number of trees to grow.
#' @param mtry Numeric. Indicates the number of variables randomly sampled as candidates at each split. This value is set sqrt(p) (where p is number of variables in data) by default.
#' @param nodesize Numeric. Indicates the minimum size of terminal nodes. Default is 1.
#' @param nvar Numeric. Indicates the number of variables to show in the Gini Index plot.
#'
#' @export
#'
#' @return A `list` with results including plots and tables.
#' @references A. Liaw and M. Wiener (2002). Classification and Regression by randomForest. R News 2(3), 18--22.
#' @author Pol Castellano-Escuder
#'
#' @importFrom magrittr %>%
#' 
#' @examples 
#' data("st000336")
#' 
#' st000336 %>% 
#'   PomaImpute() %>%
#'   PomaRandForest()
PomaRandForest <- function(data,
                           ntest = NULL,
                           ntree = 500,
                           mtry = floor(sqrt(ncol(t(SummarizedExperiment::assay(data))))),
                           nodesize = 1,
                           nvar = 20) {

  if (!is(data, "SummarizedExperiment")){
    stop("data is not a SummarizedExperiment object. \nSee POMA::PomaCreateObject or SummarizedExperiment::SummarizedExperiment")
  }
  if (!is.null(ntest)){
    if (ntest > 50 | ntest < 5) {
      stop("Incorrect value for ntest argument. It must be a number between 5 and 50")
    }
  }
  if (!is.factor(SummarizedExperiment::colData(data)[,1])) {
    stop("No factor variables to be used as dependent variable in metadata file")
  }
  
  rf_data <- data.frame(y = as.factor(as.character(SummarizedExperiment::colData(data)[,1])), 
                        t(SummarizedExperiment::assay(data))) 

  if (!is.null(ntest)) {
    
    # TRAIN AND TEST
    n <- nrow(rf_data)
    
    repeat {
      
      idx_test <- sample(1:n, (ntest/100)*n, replace = FALSE)
      
      test <- rf_data[idx_test ,]
      test_x <- as.matrix(test[,-1])
      test_y <- as.factor(test[,1])
      
      train <- rf_data[-idx_test ,]
      train_x <- as.matrix(train[,-1])
      train_y <- as.factor(train[,1])
      
      if(length(levels(as.factor(train_y))) == length(levels(as.factor(SummarizedExperiment::colData(data)[,1]))) & 
         length(levels(as.factor(test_y))) == length(levels(as.factor(SummarizedExperiment::colData(data)[,1])))) {
        break
      }
    }
    
    RF_model <- randomForest::randomForest(x = train_x,
                                           y = train_y,
                                           xtest = test_x,
                                           ytest = test_y,
                                           ntree = ntree,
                                           mtry = mtry,
                                           nodesize = nodesize)
    
  } else {
    
    RF_model <- randomForest::randomForest(y ~ ., 
                                           data = rf_data,
                                           ntree = ntree,
                                           mtry = mtry,
                                           nodesize = nodesize)
    
  }

  ntrees <- c(1:RF_model$ntree)
  error <- RF_model$err.rate

  forest_data <- data.frame(ntrees, error) %>% 
    dplyr::as_tibble()

  error_tree <- ggplot2::ggplot(forest_data, ggplot2::aes(ntrees, OOB)) +
    ggplot2::geom_line() +
    ggplot2::labs(x = "Number of Trees",
                  y = "Out-Of-Bag Error Rate") +
    theme_poma()

  importancia_pred <- randomForest::importance(RF_model, scale = TRUE) %>% 
    as.data.frame() %>% 
    dplyr::arrange(dplyr::desc(MeanDecreaseGini)) %>% 
    dplyr::slice(1:nvar) %>% 
    tibble::rownames_to_column("feature") %>% 
    dplyr::as_tibble()
  
  gini_plot <- ggplot2::ggplot(importancia_pred, ggplot2::aes(x = reorder(feature, MeanDecreaseGini),
                                                              y = MeanDecreaseGini,
                                                              fill = MeanDecreaseGini)) +
    ggplot2::labs(x = NULL,
                  y = "Mean Decrease Gini (MDG)",
                  fill = "MDG Index") +
    ggplot2::geom_col() +
    ggplot2::coord_flip() +
    theme_poma() +
    scale_fill_poma_c()

  if (!is.null(ntest)){
    cm <- caret::confusionMatrix(as.factor(RF_model$test$predicted), as.factor(test_y))
    
    return(list(MeanDecreaseGini = importancia_pred,
                MeanDecreaseGini_plot = gini_plot,
                oob_error = forest_data,
                error_tree = error_tree,
                model = RF_model,
                confusionMatrix = cm,
                train_x = train_x,
                train_y = train_y,
                test_x = test_x,
                test_y = test_y)
           )
    
  } else {
    
    return(list(MeanDecreaseGini = importancia_pred,
                MeanDecreaseGini_plot = gini_plot,
                oob_error = forest_data,
                error_tree = error_tree,
                model = RF_model)
           )
  }
}

