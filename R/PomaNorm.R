
#' Sample Quantile Normalization
#'
#' Compute quantile normalization.
#' 
#' @param data A data matrix (samples in rows).
quantile_norm <- function(data) {
  df_rank <- apply(data, 1, dplyr::dense_rank)
  df_sorted <- data.frame(apply(data, 1, sort))
  df_mean <- apply(df_sorted, 2, mean)

  index_to_mean <- function(my_index, my_mean){
    return(my_mean[my_index])
  }

  normalized_data <- apply(df_rank, 1, index_to_mean, my_mean = df_mean)
  rownames(normalized_data) <- rownames(data)
  colnames(normalized_data) <- colnames(data)
  
  return(normalized_data)
}

#' Sample Sum Normalization
#'
#' Compute sum normalization. Final unit is a percentage.
#' 
#' @param data A data matrix (samples in rows).
sum_norm <- function(data) {
  row_sums <- apply(data, 1, sum, na.rm = TRUE)
  normalized_data <- t(t(data) / row_sums) * 100
  
  return(normalized_data)
}

#' Box-Cox Transformation
#'
#' Compute Box-Cox normalization.
#' 
#' @param data A single variable.
box_cox_transformation <- function(data) {
  # Estimate optimal lambda using cross-validation
  lambda <- MASS::boxcox(data ~ 1, plotit = FALSE)
  lambda <- lambda$x[which.max(lambda$y)]
  # Perform Box-Cox transformation with the optimal lambda
  transformed <- (data^lambda - 1) / lambda
  return(transformed)
}

#' Normalize Data
#'
#' @description `PomaNorm` performs data normalization using various normalization methods.
#'
#' @param data A `SummarizedExperiment` object.
#' @param sample_norm Character. Sample normalization method. Options include "none" (default), "sum", or "quantile".
#' @param method Character. The normalization method to use. Options include "none" (no normalization), "auto_scaling" (autoscaling normalization, i.e., Z-score normalization), "level_scaling" (level scaling normalization), "log_scaling" (log scaling normalization), "log_transform" (log transformation normalization), "vast_scaling" (vast scaling normalization), "log_pareto" (log Pareto scaling normalization), "min_max" (min-max normalization), and "box_cox" (Box-Cox transformation).
#'
#' @export
#'
#' @return A `SummarizedExperiment` object with normalized data.
#' @references Van den Berg, R. A., Hoefsloot, H. C., Westerhuis, J. A., Smilde, A. K., & van der Werf, M. J. (2006). Centering, scaling, and transformations: improving the biological information content of metabolomics data. BMC genomics, 7(1), 142.
#' @author Pol Castellano-Escuder
#' 
#' @examples 
#' data("st000284")
#' 
#' PomaNorm(st000284, method = "log_pareto")
PomaNorm <- function(data,
                     sample_norm = "none",
                     method = "log_pareto"){

  if(!is(data, "SummarizedExperiment")){
    stop("data is not a SummarizedExperiment object. \nSee POMA::PomaCreateObject or SummarizedExperiment::SummarizedExperiment")
  }
  if (!(method %in% c("none", "auto_scaling", "level_scaling", "log_scaling", "log_transform",
                      "vast_scaling", "log_pareto", "min_max", "box_cox"))) {
    stop("Incorrect value for method argument")
  }
  if (missing(method)) {
    message("method argument is empty. log Pareto will be used")
  }

  to_norm <- t(SummarizedExperiment::assay(data)) %>% 
    as.data.frame()

  if (sum(is.na(to_norm)) != 0){
    stop("Missing values not allowed")
  }
  
  # remove features with only zeros
  to_norm <- data.frame(to_norm[, apply(to_norm, 2, function(x) !all(x==0))])

  # remove features with no variance
  to_norm <- data.frame(to_norm[, !apply(to_norm, 2, var) == 0])

  if (ncol(to_norm) == 1) {
    colnames(to_norm) <- t(SummarizedExperiment::assay(data)) %>% 
      as.data.frame() %>% 
      colnames()
  }
  
  # sample normalization
  if (sample_norm != "none") {
    if (sample_norm == "sum") {
      to_norm <- sum_norm(to_norm)
    }
    else if (sample_norm == "quantile") {
      to_norm <- quantile_norm(to_norm)
    }
  }
  
  # feature normalization
  if (method == "none"){
    normalized <- to_norm
  }

  else if (method == "auto_scaling"){
    normalized <- apply(to_norm, 2, function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE))
  }

  else if (method == "level_scaling"){
    normalized <- apply(to_norm, 2, function(x) (x - mean(x, na.rm = TRUE)) / mean(x, na.rm = TRUE))
  }

  else if (method == "log_scaling"){
    normalized <- apply(to_norm, 2, function(x) (log10(x + 1) - mean(log10(x + 1), na.rm = TRUE)) / sd(log10(x + 1), na.rm = TRUE))
  }

  else if (method == "log_transform"){
    normalized <- apply(to_norm, 2, function(x) (log10(x + 1)))
  }

  else if (method == "vast_scaling"){
    normalized <- apply(to_norm, 2, function(x) ((x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE))*(mean(x, na.rm = TRUE) / sd(x, na.rm = TRUE)))
  }

  else if (method == "log_pareto"){
    normalized <- apply(to_norm, 2, function(x) (log10(x + 1) - mean(log10(x + 1), na.rm = TRUE)) / sqrt(sd(log10(x + 1), na.rm = TRUE)))
  }
  
  else if (method == "min_max") {
    normalized <- apply(to_norm, 2, function(x) (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE)))
  }
  
  else if (method == "box_cox") {
    normalized <- apply(to_norm, 2, function(x) box_cox_transformation(x))
  }
  
  # create object
  if (ncol(SummarizedExperiment::colData(data)) != 0) {
    data <- SummarizedExperiment::SummarizedExperiment(assays = t(normalized), colData = SummarizedExperiment::colData(data))
  } else {
    data <- SummarizedExperiment::SummarizedExperiment(assays = t(normalized))
  }
  
  if (validObject(data))
    return(data)
}

