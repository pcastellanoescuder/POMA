
#' Sample Quantile Normalization
#'
#' Compute quantile normalization.
#' 
#' @param data A data matrix.
quantile_norm <- function(data) {
  df_rank <- apply(data, 1, dplyr::dense_rank)
  df_sorted <- data.frame(apply(data, 1, sort))
  df_mean <- apply(df_sorted, 2, mean)
  
  index_to_mean <- function(my_index, my_mean){
    return(my_mean[my_index])
  }
  
  df_final <- apply(df_rank, 1, index_to_mean, my_mean = df_mean)
  rownames(df_final) <- rownames(data)
  return(df_final)
}

#' Sample Sum Normalization
#'
#' Compute sum normalization.
#' 
#' @param data A data matrix.
sum_norm <- function(data) {
  row_sums <- apply(data, 1, sum, na.rm = TRUE)
  normalized_data <- t(t(data) / row_sums)
  
  return(normalized_data)
}

#' Collection of Normalization Methods for Mass Spectrometry Data
#'
#' @description PomaNorm() offers different methods to normalize MS data. This function contains both centering and scaling functions to normalize the data.
#'
#' @param data A SummarizedExperiment object.
#' @param method Normalization method. Options are: "none", "auto_scaling", "level_scaling", "log_scaling", "log_transformation", "vast_scaling" and "log_pareto".
#' @param sample_norm Logical. Sample sum normalization.
#' @param round Numeric. Number of decimal places (Default is 3).
#'
#' @export
#'
#' @return A SummarizedExperiment object with normalized data.
#' @references Van den Berg, R. A., Hoefsloot, H. C., Westerhuis, J. A., Smilde, A. K., & van der Werf, M. J. (2006). Centering, scaling, and transformations: improving the biological information content of metabolomics data. BMC genomics, 7(1), 142.
#' @author Pol Castellano-Escuder
#' 
#' @examples 
#' data("st000284")
#' 
#' PomaNorm(st000284, method = "log_pareto")
PomaNorm <- function(data,
                     method = "log_pareto",
                     sample_norm = TRUE,
                     round = 3){

  if (missing(data)) {
    stop("data argument is empty!")
  }
  if(!is(data, "SummarizedExperiment")){
    stop("data is not a SummarizedExperiment object. \nSee POMA::PomaSummarizedExperiment or SummarizedExperiment::SummarizedExperiment")
  }
  if (missing(method)) {
    message("method argument is empty! log_pareto will be used")
  }
  if (!(method %in% c("none", "auto_scaling", "level_scaling", "log_scaling",
                      "log_transformation", "vast_scaling", "log_pareto"))) {
    stop("Incorrect value for method argument!")
  }

  to_norm_data <- t(SummarizedExperiment::assay(data))

  # remove columns that only have zeros
  to_norm_data <- to_norm_data[, apply(to_norm_data, 2, function(x) !all(x==0))]

  # remove columns with var = 0
  to_norm_data <- to_norm_data[, !apply(to_norm_data, 2, var) == 0]

  # Sample normalization
  if (sample_norm) {
    to_norm_data <- sum_norm(to_norm_data)
  }
  
  if (method == "none"){
    normalized <- round(to_norm_data, round)
  }

  else if (method == "auto_scaling"){
    normalized <- round(apply(to_norm_data, 2, function(x) (x-mean(x,na.rm=TRUE))/sd(x,na.rm=TRUE)), round)
  }

  else if (method == "level_scaling"){
    normalized <- round(apply(to_norm_data, 2, function(x) (x-mean(x,na.rm=TRUE))/mean(x,na.rm=TRUE)), round)
  }

  else if (method == "log_scaling"){
    normalized <- round(apply(to_norm_data, 2, function(x) (log10(x+1)-mean(log10(x+1),na.rm=TRUE))/sd(log10(x+1),na.rm=TRUE)), round)
  }

  else if (method == "log_transformation"){
    normalized <- round(apply(to_norm_data, 2, function(x) (log10(x+1))), round)
  }

  else if (method == "vast_scaling"){
    normalized <- round(apply(to_norm_data, 2, function(x) ((x-mean(x,na.rm=TRUE))/sd(x,na.rm=TRUE))*(mean(x,na.rm=TRUE)/sd(x,na.rm=TRUE))), round)
  }

  else if (method == "log_pareto"){
    normalized <- round(apply(to_norm_data, 2, function(x) (log10(x+1)-mean(log10(x+1),na.rm=TRUE))/sqrt(sd(log10(x+1),na.rm=TRUE))), round)
  }
  
  target <- SummarizedExperiment::colData(data) %>% 
    as.data.frame() %>%
    tibble::rownames_to_column()
  dataNormalized <- PomaSummarizedExperiment(features = normalized, target = target)
  
  if (validObject(dataNormalized))
    return(dataNormalized)

}

