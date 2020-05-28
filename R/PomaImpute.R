
#' Collettion of Imputation Methods for Mass Spectrometry Data
#'
#' @description PomaImpute() offers different methods to impute missing values in MS data.
#'
#' @param data A MSnSet object. First `pData` column must be the subject group/type.
#' @param ZerosAsNA Logical that indicates if the zeros in the data are missing values. Default is FALSE.
#' @param RemoveNA Logical that indicates if those features with more than selected cutoff missing values in each group have to be removed. Default is TRUE.
#' @param cutoff Numeric that indicates the percentage of missing values allowed in each group. If one of the groups have less missing values than selected cutoff value, these feature will not be removed.
#' @param method Imputation method. Options are: "none", "half_min", "median", "mean", "min", "knn" and "rf". If "none", all missing values will be replaced by zero.
#'
#' @export
#'
#' @return A MSnSet object with cleaned data.
#' @references Armitage, E. G., Godzien, J., Alonso‐Herranz, V., López‐Gonzálvez, Á., & Barbas, C. (2015). Missing value imputation strategies for metabolomics data. Electrophoresis, 36(24), 3050-3060.
#' @author Pol Castellano-Escuder
#'
#' @importFrom impute impute.knn
#' @importFrom dplyr select
#' @importFrom tibble rownames_to_column
#' @importFrom magrittr %>%
#' @importFrom randomForest rfImpute
#' @importFrom crayon red
#' @importFrom clisymbols symbol
#' @importFrom Biobase varLabels pData exprs featureNames
PomaImpute <- function(data,
                       ZerosAsNA = FALSE,
                       RemoveNA = TRUE,
                       cutoff = 20,
                       method = "knn"){

  if (missing(data)) {
    stop(crayon::red(clisymbols::symbol$cross, "data argument is empty!"))
  }
  if(!(class(data) == "MSnSet")){
    stop(paste0(crayon::red(clisymbols::symbol$cross, "data is not a MSnSet object."), 
                " \nSee POMA::PomaMSnSetClass or MSnbase::MSnSet"))
  }
  if (!(method %in% c("none", "half_min", "median", "mean", "min", "knn", "rf"))) {
    stop(crayon::red(clisymbols::symbol$cross, "Incorrect value for method argument!"))
  }
  if (missing(method)) {
    warning("method argument is empty! KNN will be used")
  }

  Biobase::varLabels(data)[1] <- "Group"
  samples_groups <- Biobase::pData(data)$Group
  to_imp_data <- t(Biobase::exprs(data))
  
  ##
  
  if (isTRUE(ZerosAsNA)){
    to_imp_data[to_imp_data == 0] <- NA
    to_imp_data <- data.frame(cbind(Group = samples_groups, to_imp_data))

  } else {
    to_imp_data <- data.frame(cbind(Group = samples_groups, to_imp_data))
  }

  ##
  
  percent_na <- sum(is.na(to_imp_data))
  if (percent_na == 0) {
    stop(crayon::red(clisymbols::symbol$cross, "No missing values detected in your data"))
  }
  
  ##
  
  if (isTRUE(RemoveNA)){
    count_NA <- aggregate(. ~ Group, data = to_imp_data,
                          function(x) {100*(sum(is.na(x))/(sum(is.na(x))+sum(!is.na(x))))},
                          na.action = NULL)
    count_NA <- count_NA %>% dplyr::select(-Group)
    colnames(count_NA) <- Biobase::featureNames(data)
    supress <- as.data.frame(lapply(count_NA, function(x) all(x > cutoff)))
    supress <- unlist(supress)
    depurdata <- to_imp_data[, 2:ncol(to_imp_data)][!supress]
    depurdata <- sapply(depurdata, function(x) as.numeric(as.character(x)))

  } else {

    depurdata <- to_imp_data[, 2:ncol(to_imp_data)]
    depurdata <- sapply(depurdata, function(x) as.numeric(as.character(x)))

  }

  ##
  
  if (method == "none"){
    depurdata[is.na(depurdata)] <- 0
  }

  else if (method == "half_min"){
    depurdata <- apply(depurdata, 2, function(x) {
      if(is.numeric(x)) ifelse(is.na(x), min(x, na.rm = T)/2, x) else x})
  }

  else if (method == "median"){
    depurdata <- apply(depurdata, 2, function(x) {
      if(is.numeric(x)) ifelse(is.na(x), median(x,na.rm=T),x) else x})
  }

  else if (method == "mean"){
    depurdata <- apply(depurdata, 2, function(x) {
      if(is.numeric(x)) ifelse(is.na(x), mean(x,na.rm=T),x) else x})
  }

  else if (method == "min"){
    depurdata <- apply(depurdata, 2, function(x) {
      if(is.numeric(x)) ifelse(is.na(x), min(x,na.rm=T),x) else x})
  }

  else if (method == "knn"){
    depurdata <- t(depurdata)
    datai <- impute::impute.knn(depurdata)
    depurdata <- t(datai$data)
  }
  
  else if (method == "rf"){
    depurdata <- data.frame(group = samples_groups, depurdata)
    depurdata <- randomForest::rfImpute(group ~ ., depurdata)
    depurdata <- depurdata %>% dplyr::select(-group)
  }

  ##
  
  target <- pData(data) %>% rownames_to_column() %>% as.data.frame()
  dataImputed <- PomaMSnSetClass(features = depurdata, target = target)
  
  dataImputed@processingData@processing <-
    c(data@processingData@processing,
      paste("Imputed (", method ,"): ", date(), sep = ""))
  dataImputed@processingData@cleaned <- TRUE
  dataImputed@experimentData <- data@experimentData
  dataImputed@qual <- data@qual
    
  if (validObject(dataImputed))
    return(dataImputed)

  }

