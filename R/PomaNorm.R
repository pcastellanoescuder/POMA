
#' Collection of Normalization Methods for Mass Spectrometry Data
#'
#' @description PomaNorm() offers different methods to normalize MS data. This function contains both centering and scaling functions to normalize the data.
#'
#' @param data A MSnSet object. First `pData` column must be the subject group/type.
#' @param method Normalization method. Options are c("none", "auto_scaling", "level_scaling", "log_scaling", "log_transformation", "vast_scaling","log_pareto").
#' @param round Numeric. Number of decimal places (Default is 3).
#'
#' @export
#'
#' @return A MSnSet object with normalized data.
#' @references van den Berg, R. A., Hoefsloot, H. C., Westerhuis, J. A., Smilde, A. K., & van der Werf, M. J. (2006). Centering, scaling, and transformations: improving the biological information content of metabolomics data. BMC genomics, 7(1), 142.
#' @author Pol Castellano-Escuder
#'
#' @importFrom crayon red
#' @importFrom clisymbols symbol
#' @importFrom Biobase varLabels pData exprs
PomaNorm <- function(data,
                     method = c("none", "auto_scaling", "level_scaling", "log_scaling",
                                "log_transformation", "vast_scaling","log_pareto"),
                     round = 3){

  if (missing(method)) {
    stop(crayon::red(clisymbols::symbol$cross, "Select a method!"))
  }
  if (!(method %in% c("none", "auto_scaling", "level_scaling", "log_scaling",
                      "log_transformation", "vast_scaling","log_pareto"))) {
    stop(crayon::red(clisymbols::symbol$cross, "Incorrect value for method argument!"))
  }

  to_norm_data <- t(Biobase::exprs(data))

  # remove columns that only have zeros
  to_norm_data <- to_norm_data[, apply(to_norm_data, 2, function(x) !all(x==0))]

  # remove columns with var = 0
  to_norm_data <- to_norm_data[, !apply(to_norm_data, 2, var) == 0]

  if (method == "none"){
    normalized <- t(round(to_norm_data, round))
  }

  else if (method == "auto_scaling"){
    normalized <- t(round(apply(to_norm_data, 2, function(x) (x-mean(x,na.rm=T))/sd(x,na.rm=T)), round))
  }

  else if (method == "level_scaling"){
    normalized <- t(round(apply(to_norm_data, 2, function(x) (x-mean(x,na.rm=T))/mean(x,na.rm=T)), round))
  }

  else if (method == "log_scaling"){
    normalized <- t(round(apply(to_norm_data, 2, function(x) (log10(x+1)-mean(log10(x+1),na.rm=T))/sd(log10(x+1),na.rm=T)), round))
  }

  else if (method == "log_transformation"){
    normalized <- t(round(apply(to_norm_data, 2, function(x) (log10(x+1))), round))
  }

  else if (method == "vast_scaling"){
    normalized <- t(round(apply(to_norm_data, 2, function(x) ((x-mean(x,na.rm=T))/sd(x,na.rm=T))*(mean(x,na.rm=T)/sd(x,na.rm=T))), round))
  }

  else if (method == "log_pareto"){
    normalized <- t(round(apply(to_norm_data, 2, function(x) (log10(x+1)-mean(log10(x+1),na.rm=T))/sqrt(sd(log10(x+1),na.rm=T))), round))
  }

  data <- MSnbase::MSnSet(exprs = normalized, pData = Biobase::pData(data))
  data@processingData@processing <-
    c(data@processingData@processing,
      paste("Normalised (", method ,"): ",
            date(), sep = ""))
  data@processingData@normalised <- TRUE
  if (validObject(data))
    return(data)

}

