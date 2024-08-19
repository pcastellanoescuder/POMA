
#' Impute Missing Values
#'
#' @description `PomaImpute` performs missing value imputation on a dataset using various imputation methods.
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
#' # Output is a imputed SummarizedExperiment object
#' data <- POMA::st000284 # Example SummarizedExperiment object included in POMA
#' 
#' # No sample normalization
#' data %>% 
#'   PomaImpute(zeros_as_na = FALSE,
#'              remove_na = TRUE,
#'              cutoff = 20,
#'              group_by = TRUE,
#'              method = "knn")
PomaImpute <- function(data,
                       zeros_as_na = FALSE,
                       remove_na = TRUE,
                       cutoff = 20,
                       group_by = TRUE,
                       method = "knn"){

  if (!is(data, "SummarizedExperiment")){
    stop("data is not a SummarizedExperiment object. \nSee POMA::PomaCreateObject or SummarizedExperiment::SummarizedExperiment")
  }
  if (!(method %in% c("none", "half_min", "median", "mean", "min", "knn", "random_forest"))) {
    stop("Incorrect value for method argument")
  }
  # if (missing(method)) {
  #   message("method argument is empty. KNN will be used")
  # }

  n_features_raw <- length(rownames(data))
  
  to_impute <- t(SummarizedExperiment::assay(data)) %>% 
    as.data.frame()
  
  # zeros as NA
  if (zeros_as_na){
    to_impute[to_impute == 0] <- NA
  }
  
  percent_na <- sum(is.na(to_impute))
  if (percent_na == 0){
    message("No missing values detected")
    method <- "none"
  }
  
  grouping_factor <- ifelse(ncol(SummarizedExperiment::colData(data)) > 0, 
                            is.factor(SummarizedExperiment::colData(data)[,1]), FALSE)
  
  # remove NA
  if (remove_na){
    if (group_by & ncol(SummarizedExperiment::colData(data)) > 0 & grouping_factor) {
      to_impute <- data.frame(group_factor = SummarizedExperiment::colData(data)[,1], to_impute)
      
      count_na <- aggregate(. ~ group_factor, data = to_impute,
                            function(x) {100*(sum(is.na(x))/(sum(is.na(x)) + sum(!is.na(x))))},
                            na.action = NULL) %>%
        dplyr::select(-group_factor)
      
      to_impute <- to_impute %>%
        dplyr::select(-group_factor)
      
    } else {
      count_na <- apply(to_impute, 2, function(x) {100*(sum(is.na(x))/(sum(is.na(x)) + sum(!is.na(x))))})
    }
    remove <- unlist(as.data.frame(lapply(count_na, function(x) all(x > cutoff))))
    remove_names <- names(remove)[remove]
    
    to_impute <- to_impute %>% 
      dplyr::select(-dplyr::all_of(remove_names)) %>% 
      dplyr::mutate_all(~ as.numeric(as.character(.)))
  }
  
  # imputation
  if (method == "none"){
    to_impute[is.na(to_impute)] <- 0
    imputed <- to_impute
  }

  else if (method == "half_min"){
    imputed <- to_impute %>% 
      dplyr::mutate_all(~ ifelse(is.na(.), min(., na.rm = TRUE) / 2, .))
  }

  else if (method == "median"){
    imputed <- to_impute %>% 
      dplyr::mutate_all(~ ifelse(is.na(.), median(., na.rm = TRUE), .))
  }

  else if (method == "mean"){
    imputed <- to_impute %>% 
      dplyr::mutate_all(~ ifelse(is.na(.), mean(., na.rm = TRUE), .))
  }

  else if (method == "min"){
    imputed <- to_impute %>% 
      dplyr::mutate_all(~ ifelse(is.na(.), min(., na.rm = TRUE), .))
  }

  else if (method == "knn"){
    suppressWarnings({
      imputed_t <- t(to_impute)
      imputed_res <- impute::impute.knn(imputed_t)
      imputed <- t(imputed_res$data)
    })
  }
  
  else if (method == "random_forest"){
    if (ncol(SummarizedExperiment::colData(data)) == 0 | !is.factor(SummarizedExperiment::colData(data)[,1])) {
      stop("This imputation method is not compatible with the provided metadata")
    }
    imputed <- data.frame(group_factor = SummarizedExperiment::colData(data)[,1], to_impute)
    imputed <- randomForest::rfImpute(group_factor ~ ., imputed)
    imputed <- imputed %>% 
      dplyr::select(-group_factor)
  }
  
  # create object
  if (ncol(SummarizedExperiment::colData(data)) != 0) {
    data <- SummarizedExperiment::SummarizedExperiment(assays = t(imputed), colData = SummarizedExperiment::colData(data))
  } else {
    data <- SummarizedExperiment::SummarizedExperiment(assays = t(imputed))
  }
  
  n_features_imputed <- length(rownames(data))
  
  message(paste0(n_features_raw - n_features_imputed, " features removed."))
  
  if (validObject(data))
    return(data)
}

