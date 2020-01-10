
#' Classification Random Forest for Mass Spectrometry Data
#'
#' @description PomaRandForest() allows users to perform a classification Random Forest with a MS data matrix using the classical `randomForest` R package.
#'
#' @param data A MSnSet object. First `pData` column must be the suject group/type.
#' @param folds Number of observations that will be used as test dataset. For example, if folds = 3, 1/3 of dataset will be used as train dataset.
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
#' @importFrom tibble rownames_to_column
#' @importFrom crayon red
#' @importFrom clisymbols symbol
#' @importFrom Biobase varLabels pData exprs
PomaRandForest <- function(data,
                           folds = 3,
                           ntree = 500,
                           mtry = floor(sqrt(ncol(t(Biobase::exprs(data))))),
                           nodesize = 1,
                           nvar = 20){

  if (missing(data)) {
    stop(crayon::red(clisymbols::symbol$cross, "data argument is empty!"))
  }

  Biobase::varLabels(data)[1] <- "Group"
  rf_data <- data.frame(cbind(Group = Biobase::pData(data)$Group, t(Biobase::exprs(data))))
  rf_data[, 2:ncol(rf_data)] <- sapply(rf_data[, 2:ncol(rf_data)], as.numeric)

  names <- data.frame(real_names = colnames(rf_data), new_names = NA)
  names$new_names <- paste0("X", rownames(names))

  colnames(rf_data) <- names$new_names
  colnames(rf_data)[1] <- "Group"

  # training Sample with 1/folds observations
  train <- sample(1:nrow(rf_data), round(nrow(rf_data)/folds))

  RF_model <- randomForest(as.factor(Group) ~ . ,
                           data = rf_data,
                           subset = train,
                           ntree = ntree,
                           mtry = mtry,
                           nodesize = nodesize)

  ntrees <- c(1:RF_model$ntree)
  error <- RF_model$err.rate

  forest_data <- round(data.frame(ntrees, error), 4)

  error_tree <- ggplot(forest_data, aes(ntrees, OOB)) +
    geom_line() +
    labs(y = "Out-Of-Bag Error Rate") +
    theme_minimal()

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
    xlab("Variable") +
    geom_col() +
    coord_flip() +
    theme_minimal() +
    theme(legend.position = "bottom")

  importancia_pred1 <- importancia_pred1[, c(3,2)]
  importancia_pred1$MeanDecreaseGini <- round(importancia_pred1$MeanDecreaseGini, 4)
  colnames(importancia_pred1)[1] <- "Variable"

  conf_mat <- round(as.data.frame(RF_model$confusion), 4)

  return(list(importance_pred = importancia_pred1,
              error_tree = error_tree,
              gini_plot = Gini_plot,
              forest_data = forest_data,
              confusion_matrix = conf_mat))
}

