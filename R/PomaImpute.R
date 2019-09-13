
#' Different Imputation Methods for Metabolomics
#'
#' @description PomaImpute() offers different methods to impute missing values in metabolomic data.
#'
#' @param data A data frame with metabolites. First column must be the subject ID and second column must be a factor with the subject group.
#' @param ZerosAsNA Logical that indicates if the zeros in the data are missing values. Default is FALSE.
#' @param RemoveNA Logical that indicates if those metabolites with more than selected cutoff missing values in each group have to be removed. Default is TRUE.
#' @param cutoff Numeric that indicates the percentage of missing values allowed in each group. If one of the groups have less missing values than selected cutoff value, these metabolite will not be removed.
#' @param method Imputation method. Options are c("none", "half_min", "median", "mean", "min", "knn"). If "none", all missing values will be replaced by zero.
#'
#' @export
#'
#' @return A data frame with the results.
#' @references Armitage, E. G., Godzien, J., Alonso‐Herranz, V., López‐Gonzálvez, Á., & Barbas, C. (2015). Missing value imputation strategies for metabolomics data. Electrophoresis, 36(24), 3050-3060.
#' @author Pol Castellano-Escuder
PomaImpute <- function(data,
                       ZerosAsNA = FALSE,
                       RemoveNA = TRUE,
                       cutoff = 20,
                       method = c("none", "half_min", "median", "mean", "min", "knn")){

  if (!(method %in% c("none", "half_min", "median", "mean", "min", "knn"))) {
    stop(crayon::red(clisymbols::symbol$cross, "Incorrect value for method argument!"))
  }
  if (missing(method)) {
    method <- "knn"
    warning("method argument is empty! KNN will be used")
  }

  data <- as.data.frame(data)
  samples_groups <- data[, 1:2]
  to_imp_data <- data[, c(3:ncol(data))]
  colnames(samples_groups)[2] <- "Group"

  if (ZerosAsNA == TRUE){
    to_imp_data[to_imp_data == 0] <- NA
    to_imp_data <- cbind(samples_groups, to_imp_data)
  }

  else{
    to_imp_data <- cbind(samples_groups, to_imp_data)
  }

  if (RemoveNA == TRUE){
    count_NA <- aggregate(. ~ Group, data = to_imp_data[, 2:ncol(to_imp_data)],
                          function(x) {100*(sum(is.na(x))/(sum(is.na(x))+sum(!is.na(x))))},
                          na.action = NULL)
    count_NA$Group <- NULL
    supress <- as.data.frame(lapply(count_NA, function(x) all(x > cutoff)))
    supress <- unlist(supress)
    depurdata <- to_imp_data[, 3:ncol(to_imp_data)][!supress]

    if (method == "none"){
      depurdata[is.na(depurdata)] <- 0
      depurdata <- cbind(samples_groups, depurdata)
      return(depurdata)
    }

    else if (method == "half_min"){
      depurdata <- apply(depurdata, 2, function(x) "[<-"(x, !x | is.na(x),
                                                         min(x[x >= 0], na.rm = TRUE) / 2))
      depurdata <- cbind(samples_groups, depurdata)
      return(depurdata)
    }

    else if (method == "median"){
      depurdata <- apply(depurdata, 2, function(x) {
        if(is.numeric(x)) ifelse(is.na(x), median(x,na.rm=T),x) else x})
      depurdata <- cbind(samples_groups, depurdata)
      return(depurdata)
    }

    else if (method == "mean"){
      depurdata <- apply(depurdata, 2, function(x) {
        if(is.numeric(x)) ifelse(is.na(x), mean(x,na.rm=T),x) else x})
      depurdata <- cbind(samples_groups, depurdata)
      return(depurdata)
    }

    else if (method == "min"){
      depurdata <- apply(depurdata, 2, function(x) {
        if(is.numeric(x)) ifelse(is.na(x), min(x,na.rm=T),x) else x})
      depurdata <- cbind(samples_groups, depurdata)
      return(depurdata)
    }

    else if (method == "knn"){
      depurdata <- t(depurdata)
      datai <- impute::impute.knn(as.matrix(depurdata))
      depurdata <- t(datai$data)
      depurdata <- cbind(samples_groups, depurdata)
      return(depurdata)
    }
    return (depurdata)

  }

  else{

    depurdata <- to_imp_data[, 3:ncol(to_imp_data)]

    if (method == "none"){
      depurdata[is.na(depurdata)] <- 0
      depurdata <- cbind(samples_groups, depurdata)
      return(depurdata)
    }

    else if (method == "half_min"){
      depurdata <- apply(depurdata, 2, function(x) "[<-"(x, !x | is.na(x),
                                                         min(x[x >= 0], na.rm = TRUE) / 2))
      depurdata <- cbind(samples_groups, depurdata)
      return(depurdata)
    }

    else if (method == "median"){
      depurdata <- apply(depurdata, 2, function(x) {
        if(is.numeric(x)) ifelse(is.na(x), median(x,na.rm=T),x) else x})
      depurdata <- cbind(samples_groups, depurdata)
      return(depurdata)
    }

    else if (method == "mean"){
      depurdata <- apply(depurdata, 2, function(x) {
        if(is.numeric(x)) ifelse(is.na(x), mean(x,na.rm=T),x) else x})
      depurdata <- cbind(samples_groups, depurdata)
      return(depurdata)
    }

    else if (method == "min"){
      depurdata <- apply(depurdata, 2, function(x) {
        if(is.numeric(x)) ifelse(is.na(x), min(x,na.rm=T),x) else x})
      depurdata <- cbind(samples_groups, depurdata)
      return(depurdata)
    }

    else if (method == "knn"){
      depurdata <- t(depurdata)
      datai <- impute::impute.knn(as.matrix(depurdata))
      depurdata <- t(datai$data)
      depurdata <- cbind(samples_groups, depurdata)
      return(depurdata)
    }
    return (depurdata)
  }

}

