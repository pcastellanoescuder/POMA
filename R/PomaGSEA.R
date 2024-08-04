
#' Gene Set Enrichment Analysis
#'
#' @description `PomaGSEA` performs missing value imputation on a dataset using various imputation methods.
#'
#' @param data A `SummarizedExperiment` object.
#' @param zeros_as_na Logical. Indicates if the zeros in the data are missing values. Default is FALSE.
#' @param remove_na Logical. Indicates if features with a percentage of missing values over the `cutoff` parameter should be removed. Default is TRUE.
#' @param cutoff Numeric. Percentage of missing values allowed in each feature.
#' @param group_by Logical. If `metadata` file is present and its first variable is a factor, it can be used to compute missing values per group and drop them accordingly. Features will be removed only if all of the groups contain more missing values than allowed. Default is TRUE.
#' @param method Character. The imputation method to use. Options include "none" (no imputation, replace missing values by zeros), "half_min" (replace missing values with half of the minimum value), "median" (replace missing values with the median), "mean" (replace missing values with the mean), "min" (replace missing values with the minimum value), "knn" (replace missing values using k-nearest neighbors imputation), and "random_forest" (replace missing values using random forest imputation).
#'
#' @export
#'
#' @return A `SummarizedExperiment` object without missing values.
#' @references Armitage, E. G., Godzien, J., Alonso‐Herranz, V., López‐Gonzálvez, Á., & Barbas, C. (2015). Missing value imputation strategies for metabolomics data. Electrophoresis, 36(24), 3050-3060.
#' @author Pol Castellano-Escuder
#'
#' @importFrom magrittr %>%
#' 
#' @examples 
#' data("st000336")
#' 
#' PomaGSEA(st000336, method = "knn")
PomaGSEA <- function(data) {
  
  ranked_data <- data %>% 
    as.data.frame() %>% 
    dplyr::select(feature = 1, rank = 2)
  
}

