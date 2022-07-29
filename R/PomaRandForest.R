
#' Classification Random Forest
#'
#' @description PomaRandForest() allows users to perform a classification Random Forest using the classical `randomForest` R package.
#'
#' @param data A SummarizedExperiment object.
#' @param ntest Numeric indicating the percentage of observations that will be used as test set. Default is NULL (no test set).
#' @param ntree Number of trees to grow.
#' @param mtry Number of variables randomly sampled as candidates at each split. This value is set sqrt(p) (where p is number of variables in data) by default.
#' @param nodesize Minimum size of terminal nodes. By default is equal to 1.
#' @param nvar Number of variables to show in the Gini plot.
#'
#' @export
#'
#' @return A list with all results for Random Forest including plots and tables.
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
                           nvar = 20){

  if (missing(data)) {
    stop("data argument is empty!")
  }
  if(!is(data, "SummarizedExperiment")){
    stop("data is not a SummarizedExperiment object. \nSee POMA::PomaSummarizedExperiment or SummarizedExperiment::SummarizedExperiment")
  }
  if(!is.null(ntest)){
    if (ntest > 50 | ntest < 5) {
      stop("ntest must be a number between 5 and 50...")
    }
  }
  
  rf_data <- data.frame(cbind(Group = SummarizedExperiment::colData(data)[,1], 
                              t(SummarizedExperiment::assay(data))))

  names <- data.frame(real_names = colnames(rf_data), new_names = NA) %>%
    dplyr::mutate(new_names = paste0("X", rownames(.)))

  colnames(rf_data) <- names$new_names
  colnames(rf_data)[1] <- "Group"

  if(!is.null(ntest)){
    
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
         length(levels(as.factor(test_y))) == length(levels(as.factor(SummarizedExperiment::colData(data)[,1])))){
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
    
    rf_data <- rf_data %>% 
      dplyr::mutate(Group = as.factor(Group))
    
    RF_model <- randomForest::randomForest(Group ~ ., 
                                           data = rf_data,
                                           ntree = ntree,
                                           mtry = mtry,
                                           nodesize = nodesize)
    
  }

  ntrees <- c(1:RF_model$ntree)
  error <- RF_model$err.rate

  forest_data <- round(data.frame(ntrees, error), 4) %>% 
    dplyr::as_tibble()

  error_tree <- ggplot2::ggplot(forest_data, ggplot2::aes(ntrees, OOB)) +
    ggplot2::geom_line() +
    ggplot2::labs(y = "Out-Of-Bag Error Rate") +
    ggplot2::theme_bw()

  importancia_pred <- randomForest::importance(RF_model, scale = TRUE) %>% 
    as.data.frame() %>% 
    tibble::rownames_to_column("new_names") %>% 
    dplyr::inner_join(names, by = "new_names") %>% 
    dplyr::arrange(dplyr::desc(MeanDecreaseGini)) %>% 
    dplyr::slice(1:nvar) %>% 
    dplyr::select(real_names, MeanDecreaseGini) %>% 
    dplyr::mutate(MeanDecreaseGini = round(MeanDecreaseGini, 3)) %>% 
    dplyr::rename(feature = real_names) %>% 
    dplyr::as_tibble()
  
  gini_plot <- ggplot2::ggplot(importancia_pred, ggplot2::aes(x = reorder(feature, MeanDecreaseGini),
                                                              y = MeanDecreaseGini,
                                                              fill = MeanDecreaseGini)) +
    ggplot2::labs(x = NULL) +
    ggplot2::geom_col() +
    ggplot2::coord_flip() +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "none") +
    ggplot2::scale_fill_viridis_c(begin = 0, end = 0.8)

  feature_names <- names %>% 
    dplyr::slice(-1) %>% 
    dplyr::rename(feature = real_names,
                  idx = new_names) %>% 
    dplyr::as_tibble()

  if(!is.null(ntest)){
    cm <- caret::confusionMatrix(as.factor(RF_model$test$predicted), as.factor(test_y))
    
    return(list(MeanDecreaseGini = importancia_pred,
                MeanDecreaseGini_plot = gini_plot,
                oob_error = forest_data,
                error_tree = error_tree,
                feature_names = feature_names,
                model = RF_model,
                confusionMatrix = cm,
                train_x = train_x,
                train_y = train_y,
                test_x = test_x,
                test_y = test_y))
    
  } else {
    
    return(list(MeanDecreaseGini = importancia_pred,
                MeanDecreaseGini_plot = gini_plot,
                oob_error = forest_data,
                error_tree = error_tree,
                feature_names = feature_names,
                model = RF_model))
  }
  
}

