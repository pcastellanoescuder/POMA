
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

  Biobase::varLabels(st000284)[1] <- "Group"
  samples_groups <- Biobase::pData(st000284)$Group
  to_imp_data <- t(Biobase::exprs(st000284))

  if (ZerosAsNA == TRUE){
    to_imp_data[to_imp_data == 0] <- NA
    to_imp_data <- data.frame(cbind(Group = samples_groups, to_imp_data))
    # to_imp_data[,2:ncol(to_imp_data)] <- sapply(to_imp_data[,2:ncol(to_imp_data)], as.numeric)
  }

  else{
    to_imp_data <- data.frame(cbind(Group = samples_groups, to_imp_data))
    # to_imp_data[,2:ncol(to_imp_data)] <- sapply(to_imp_data[,2:ncol(to_imp_data)], as.numeric)
  }

  if (RemoveNA == TRUE){
    count_NA <- aggregate(. ~ Group, data = to_imp_data,
                          function(x) {100*(sum(is.na(x))/(sum(is.na(x))+sum(!is.na(x))))},
                          na.action = NULL)
    count_NA$Group <- NULL
    supress <- as.data.frame(lapply(count_NA, function(x) all(x > cutoff)))
    supress <- unlist(supress)
    depurdata <- to_imp_data[, 2:ncol(to_imp_data)][!supress]

    if (method == "none"){
      depurdata[is.na(depurdata)] <- 0
    }

    else if (method == "half_min"){
      depurdata <- apply(depurdata, 2, function(x) "[<-"(x, !x | is.na(x),
                                                         min(x[x >= 0], na.rm = TRUE) / 2))
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
      datai <- impute::impute.knn(as.matrix(depurdata))
      depurdata <- t(datai$data)
    }

  }

  else{

    if (method == "none"){
      to_imp_data[is.na(to_imp_data)] <- 0
      depurdata <- to_imp_data
    }

    else if (method == "half_min"){
      depurdata <- apply(to_imp_data, 2, function(x) "[<-"(x, !x | is.na(x),
                                                         min(x[x >= 0], na.rm = TRUE) / 2))
    }

    else if (method == "median"){
      depurdata <- apply(to_imp_data, 2, function(x) {
        if(is.numeric(x)) ifelse(is.na(x), median(x,na.rm=T),x) else x})
    }

    else if (method == "mean"){
      depurdata <- apply(to_imp_data, 2, function(x) {
        if(is.numeric(x)) ifelse(is.na(x), mean(x,na.rm=T),x) else x})
    }

    else if (method == "min"){
      depurdata <- apply(to_imp_data, 2, function(x) {
        if(is.numeric(x)) ifelse(is.na(x), min(x,na.rm=T),x) else x})
    }

    else if (method == "knn"){
      depurdata <- t(to_imp_data)
      datai <- impute::impute.knn(as.matrix(depurdata))
      depurdata <- t(datai$data)
    }
  }

  Biobase::exprs(data) <- depurdata
  # data@processingData@processing <-
  #   c(data@processingData@processing,
  #     paste("Imputed (", method ,"): ",
  #           date(), sep = ""))
  # data@processingData@cleaned <- TRUE
  # if (validObject(data))
    return(data)

}

