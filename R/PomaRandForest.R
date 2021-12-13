
#' Classification Random Forest for Mass Spectrometry Data
#'
#' @description PomaRandForest() allows users to perform a classification Random Forest with a MS data matrix using the classical `randomForest` R package.
#'
#' @param data A SummarizedExperiment object. First `colData` column must be the subject group/type.
#' @param ntest Numeric indicating the percentage of observations that will be used as test set. Default is 20% of observations.
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
#' @importFrom randomForest randomForest importance
#' @import ggplot2
#' @importFrom dplyr mutate as_tibble inner_join arrange desc slice select rename
#' @importFrom magrittr %>%
#' @importFrom tibble rownames_to_column
#' @importFrom SummarizedExperiment assay colData
#' 
#' @examples 
#' data("st000336")
#' 
#' st000336 %>% 
#'   PomaImpute() %>%
#'   PomaRandForest()
PomaRandForest <- function(data,
                           ntest = 20,
                           ntree = 500,
                           mtry = floor(sqrt(ncol(t(SummarizedExperiment::assay(data))))),
                           nodesize = 1,
                           nvar = 20){

  if (missing(data)) {
    stop("data argument is empty!")
  }
  if(!is(data[1], "SummarizedExperiment")){
    stop("data is not a SummarizedExperiment object. \nSee POMA::PomaSummarizedExperiment or SummarizedExperiment::SummarizedExperiment")
  }
  if (ntest > 50 | ntest < 5) {
    stop("ntest must be a number between 5 and 50...")
  }
  
  rf_data <- data.frame(cbind(Group = SummarizedExperiment::colData(data)[,1], 
                              t(SummarizedExperiment::assay(data))))

  names <- data.frame(real_names = colnames(rf_data), new_names = NA) %>%
    mutate(new_names = paste0("X", rownames(.)))

  colnames(rf_data) <- names$new_names
  colnames(rf_data)[1] <- "Group"

  # TRAIN AND TEST
  n <- nrow(rf_data)
  
  repeat{
    
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

  RF_model <- randomForest(x = train_x,
                           y = train_y,
                           xtest = test_x,
                           ytest = test_y,
                           ntree = ntree,
                           mtry = mtry,
                           nodesize = nodesize)

  ntrees <- c(1:RF_model$ntree)
  error <- RF_model$err.rate

  forest_data <- round(data.frame(ntrees, error), 4) %>% 
    as_tibble()

  error_tree <- ggplot(forest_data, aes(ntrees, OOB)) +
    geom_line() +
    labs(y = "Out-Of-Bag Error Rate") +
    theme_bw()

  ####

  importancia_pred <- randomForest::importance(RF_model, scale = TRUE) %>% 
    as.data.frame() %>% 
    rownames_to_column("new_names") %>% 
    inner_join(names, by = "new_names") %>% 
    arrange(desc(MeanDecreaseGini)) %>% 
    dplyr::slice(1:nvar) %>% 
    dplyr::select(real_names, MeanDecreaseGini) %>% 
    dplyr::mutate(MeanDecreaseGini = round(MeanDecreaseGini, 3)) %>% 
    dplyr::rename(feature = real_names) %>% 
    dplyr::as_tibble()

  ##
  
  gini_plot <- ggplot(importancia_pred, aes(x = reorder(feature, MeanDecreaseGini),
                                            y = MeanDecreaseGini,
                                            fill = MeanDecreaseGini)) +
    xlab("") +
    geom_col() +
    coord_flip() +
    theme_bw() +
    theme(legend.position = "none")
  
  ##

  feature_names <- names %>% 
    dplyr::slice(-1) %>% 
    dplyr::rename(feature = real_names,
                  idx = new_names) %>% 
    dplyr::as_tibble()
    
  ##
  
  conf_mat <- round(as.data.frame(RF_model$test$confusion), 4)

  return(list(MeanDecreaseGini = importancia_pred,
              MeanDecreaseGini_plot = gini_plot,
              oob_error = forest_data,
              error_tree = error_tree,
              confusion_matrix = conf_mat,
              feature_names = feature_names,
              model = RF_model,
              train_x = train_x,
              train_y = train_y,
              test_x = test_x,
              test_y = test_y))
}

