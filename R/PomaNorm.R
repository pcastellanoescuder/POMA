
#' Collection of Normalization Methods for Mass Spectrometry Data
#'
#' @description PomaNorm() offers different methods to normalize MS data. This function contains both centering and scaling functions to normalize the data.
#'
#' @param data A MSnSet object. First `pData` column must be the subject group/type.
#' @param method Normalization method. Options are: "none", "auto_scaling", "level_scaling", "log_scaling", "log_transformation", "vast_scaling" and "log_pareto".
#' @param round Numeric. Number of decimal places (Default is 3).
#'
#' @export
#'
#' @return A MSnSet object with normalized data.
#' @references van den Berg, R. A., Hoefsloot, H. C., Westerhuis, J. A., Smilde, A. K., & van der Werf, M. J. (2006). Centering, scaling, and transformations: improving the biological information content of metabolomics data. BMC genomics, 7(1), 142.
#' @author Pol Castellano-Escuder
#'
#' @importFrom crayon red
#' @importFrom tibble rownames_to_column
#' @importFrom clisymbols symbol
#' @importFrom Biobase varLabels pData exprs
#' 
#' @examples 
#' data("st000284")
#' 
#' PomaNorm(st000284, method = "log_pareto")
PomaNorm <- function(data,
                     method = "log_pareto",
                     round = 3){

  if (missing(data)) {
    stop(crayon::red(clisymbols::symbol$cross, "data argument is empty!"))
  }
  if(!is(data[1], "MSnSet")){
    stop(paste0(crayon::red(clisymbols::symbol$cross, "data is not a MSnSet object."), 
                " \nSee POMA::PomaMSnSetClass or MSnbase::MSnSet"))
  }
  if (missing(method)) {
    warning("method argument is empty! log_pareto will be used")
  }
  if (!(method %in% c("none", "auto_scaling", "level_scaling", "log_scaling",
                      "log_transformation", "vast_scaling", "log_pareto"))) {
    stop(crayon::red(clisymbols::symbol$cross, "Incorrect value for method argument!"))
  }

  to_norm_data <- t(Biobase::exprs(data))

  # remove columns that only have zeros
  to_norm_data <- to_norm_data[, apply(to_norm_data, 2, function(x) !all(x==0))]

  # remove columns with var = 0
  to_norm_data <- to_norm_data[, !apply(to_norm_data, 2, var) == 0]

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

  ##
  
  target <- pData(data) %>% rownames_to_column() %>% as.data.frame()
  dataNormalized <- PomaMSnSetClass(features = normalized, target = target)
  
  dataNormalized@processingData@processing <-
    c(data@processingData@processing,
      paste("Normalised (", method ,"): ", date(), sep = ""))
  dataNormalized@processingData@normalised <- TRUE
  dataNormalized@processingData@cleaned <- data@processingData@cleaned
  dataNormalized@experimentData <- data@experimentData
  dataNormalized@qual <- data@qual
  
  if (validObject(dataNormalized))
    return(dataNormalized)

}

