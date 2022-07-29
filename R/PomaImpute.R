
#' Collection of Imputation Methods for Mass Spectrometry Data
#'
#' @description PomaImpute() offers different methods to impute missing values in MS data.
#'
#' @param data A SummarizedExperiment object.
#' @param ZerosAsNA Logical that indicates if the zeros in the data are missing values. Default is FALSE.
#' @param RemoveNA Logical that indicates if those features with more than selected cutoff missing values in each group have to be removed. Default is TRUE.
#' @param cutoff Numeric that indicates the percentage of missing values allowed in each group. If one of the groups have less missing values than selected cutoff value, these feature will not be removed.
#' @param method Imputation method. Options are: "none", "half_min", "median", "mean", "min", "knn" and "rf". If "none", all missing values will be replaced by zero.
#'
#' @export
#'
#' @return A SummarizedExperiment object with cleaned data.
#' @references Armitage, E. G., Godzien, J., Alonso‐Herranz, V., López‐Gonzálvez, Á., & Barbas, C. (2015). Missing value imputation strategies for metabolomics data. Electrophoresis, 36(24), 3050-3060.
#' @author Pol Castellano-Escuder
#'
#' @importFrom magrittr %>%
#' 
#' @examples 
#' data("st000336")
#' 
#' PomaImpute(st000336, method = "knn")
PomaImpute <- function(data,
                       ZerosAsNA = FALSE,
                       RemoveNA = TRUE,
                       cutoff = 20,
                       method = "knn"){

  if (missing(data)) {
    stop("data argument is empty!")
  }
  if(!is(data, "SummarizedExperiment")){
    stop("data is not a SummarizedExperiment object. \nSee POMA::PomaSummarizedExperiment or SummarizedExperiment::SummarizedExperiment")
  }
  if (!(method %in% c("none", "half_min", "median", "mean", "min", "knn", "rf"))) {
    stop("Incorrect value for method argument!")
  }
  if (missing(method)) {
    message("method argument is empty! KNN will be used")
  }

  samples_groups <- SummarizedExperiment::colData(data)[,1]
  to_imp_data <- t(SummarizedExperiment::assay(data))
  
  ##
  
  if (ZerosAsNA){
    to_imp_data[to_imp_data == 0] <- NA
    to_imp_data <- data.frame(cbind(Group = samples_groups, to_imp_data))
    colnames(to_imp_data)[2:ncol(to_imp_data)] <- data@NAMES

  } else {
    to_imp_data <- data.frame(cbind(Group = samples_groups, to_imp_data))
    colnames(to_imp_data)[2:ncol(to_imp_data)] <- data@NAMES
  }

  ##
  
  percent_na <- sum(is.na(to_imp_data))
  if (percent_na == 0){
    message("No missing values detected in your data")
    if(method == "rf") {
      method <- "none"
    }
  }
  
  ##
  
  if (isTRUE(RemoveNA)){
    count_NA <- aggregate(. ~ Group, data = to_imp_data,
                          function(x) {100*(sum(is.na(x))/(sum(is.na(x))+sum(!is.na(x))))},
                          na.action = NULL)
    count_NA <- count_NA %>% 
      dplyr::select(-Group)
    correct_names <- names(count_NA)
    supress <- unlist(as.data.frame(lapply(count_NA, function(x) all(x > cutoff))))
    names(supress) <- correct_names
    correct_names <- names(supress[supress == "FALSE"])
    depurdata <- to_imp_data[, 2:ncol(to_imp_data)][!supress]
    depurdata <- sapply(depurdata, function(x) as.numeric(as.character(x)))

  } else {
    depurdata <- to_imp_data[, 2:ncol(to_imp_data)]
    depurdata <- sapply(depurdata, function(x) as.numeric(as.character(x)))
    correct_names <- data@NAMES
  }

  ##
  
  if (method == "none"){
    depurdata[is.na(depurdata)] <- 0
  }

  else if (method == "half_min"){
    depurdata <- apply(depurdata, 2, function(x) {
      if(is.numeric(x)) ifelse(is.na(x), min(x, na.rm = TRUE)/2, x) else x})
  }

  else if (method == "median"){
    depurdata <- apply(depurdata, 2, function(x) {
      if(is.numeric(x)) ifelse(is.na(x), median(x, na.rm = TRUE), x) else x})
  }

  else if (method == "mean"){
    depurdata <- apply(depurdata, 2, function(x) {
      if(is.numeric(x)) ifelse(is.na(x), mean(x, na.rm = TRUE), x) else x})
  }

  else if (method == "min"){
    depurdata <- apply(depurdata, 2, function(x) {
      if(is.numeric(x)) ifelse(is.na(x), min(x, na.rm = TRUE), x) else x})
  }

  else if (method == "knn"){
    depurdata <- t(depurdata)
    datai <- impute::impute.knn(depurdata)
    depurdata <- t(datai$data)
  }
  
  else if (method == "rf"){
    depurdata <- data.frame(group = samples_groups, depurdata)
    depurdata <- randomForest::rfImpute(group ~ ., depurdata)
    depurdata <- depurdata %>% 
      dplyr::select(-group)
  }
  
  colnames(depurdata) <- correct_names
  
  target <- SummarizedExperiment::colData(data) %>% 
    as.data.frame() %>%
    tibble::rownames_to_column()
  dataImputed <- PomaSummarizedExperiment(features = depurdata, target = target)
    
  if (validObject(dataImputed))
    return(dataImputed)

}

