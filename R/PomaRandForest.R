
#' Classification Random Forest for Mass Spectrometry Data
#'
#' @description PomaRandForest() allows users to perform a classification Random Forest with a MS data matrix using the classical `randomForest` R package.
#'
#' @param data A MSnSet object. First `pData` column must be the subject group/type.
#' @param ntest Numeric indicating the percentage of observations that will be used as test set. Default is 20% of observations.
#' @param ntree Number of trees to grow.
#' @param mtry Number of variables randomly sampled as candidates at each split. This value is set sqrt(p) (where p is number of variables in data) by default.
#' @param nodesize Minimum size of terminal nodes. By default is equal to 1.
#' @param nvar Number of variables to show in the Gini plot.
#'
#' @export
#'
#' @return A list with all results for Random Forest including plots and data frames.
#' @references A. Liaw and M. Wiener (2002). Classification and Regression by randomForest. R News 2(3), 18--22.
#' @author Pol Castellano-Escuder
#'
#' @importFrom randomForest randomForest importance
#' @import ggplot2
#' @importFrom dplyr mutate
#' @importFrom magrittr %>%
#' @importFrom tibble rownames_to_column
#' @importFrom crayon red
#' @importFrom clisymbols symbol
#' @importFrom Biobase varLabels pData exprs
PomaRandForest <- function(data,
                           ntest = 20,
                           ntree = 500,
                           mtry = floor(sqrt(ncol(t(Biobase::exprs(data))))),
                           nodesize = 1,
                           nvar = 20){

  if (missing(data)) {
    stop(crayon::red(clisymbols::symbol$cross, "data argument is empty!"))
  }
  if(!(class(data) == "MSnSet")){
    stop(paste0(crayon::red(clisymbols::symbol$cross, "data is not a MSnSet object."), 
                " \nSee POMA::PomaMSnSetClass or MSnbase::MSnSet"))
  }
  if (ntest > 50 | ntest < 5) {
    stop(crayon::red(clisymbols::symbol$cross, "ntest must be a number between 5 and 50..."))
  }
  

  Biobase::varLabels(data)[1] <- "Group"
  rf_data <- data.frame(cbind(Group = Biobase::pData(data)$Group, t(Biobase::exprs(data))))

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
    
    if(length(levels(as.factor(train_y))) == length(levels(as.factor(Biobase::pData(data)[,1]))) & 
       length(levels(as.factor(test_y))) == length(levels(as.factor(Biobase::pData(data)[,1])))){
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

  forest_data <- round(data.frame(ntrees, error), 4)

  error_tree <- ggplot(forest_data, aes(ntrees, OOB)) +
    geom_line() +
    labs(y = "Out-Of-Bag Error Rate") +
    theme_bw()

  ####

  importancia_pred <- as.data.frame(randomForest::importance(RF_model, scale = TRUE))
  importancia_pred <- rownames_to_column(importancia_pred, var = "new_names")
  importancia_pred <- merge(importancia_pred, names , by = "new_names")
  importancia_pred1 <- importancia_pred[order(importancia_pred$MeanDecreaseGini,
                                              decreasing = TRUE),]
  importancia_pred <- importancia_pred1[1:nvar ,]

  Gini_plot <- ggplot(importancia_pred, aes(x = reorder(real_names, MeanDecreaseGini),
                                            y = MeanDecreaseGini,
                                            fill = MeanDecreaseGini)) +
    xlab("") +
    geom_col() +
    coord_flip() +
    theme_bw() +
    theme(legend.position = "none")

  importancia_pred1 <- importancia_pred1[, c(3,2)]
  importancia_pred1$MeanDecreaseGini <- round(importancia_pred1$MeanDecreaseGini, 4)
  colnames(importancia_pred1)[1] <- "Variable"

  conf_mat <- round(as.data.frame(RF_model$test$confusion), 4)

  return(list(importance_pred = importancia_pred1,
              error_tree = error_tree,
              gini_plot = Gini_plot,
              forest_data = forest_data,
              confusion_matrix = conf_mat,
              model = RF_model))
}

