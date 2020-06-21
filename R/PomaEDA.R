
#' Automatic Exploratory Data Analysis PDF Report
#'
#' @description This function automatically generates a PDF report with different exploratory plots and tables from an MSnSet object.
#'
#' @param data A MSnSet object. First `pData` column must be the subject group/type.
#' @param imputation Imputation method. Options are "none", "half_min", "median", "mean", "min" and "knn" (default). If "none", all missing values will be replaced by zero.
#' @param normalization Normalization method. Options are "none", "auto_scaling", "level_scaling", "log_scaling", "log_transformation", "vast_scaling" and "log_pareto" (default).
#' @param clean_outliers Logical. If it's set to TRUE, outliers will be removed from EDA.
#' @param coeff_outliers This value corresponds to the classical 1.5 in \eqn{Q3 + 1.5*IQR} formula to detect outliers. By changing this value, the permissiveness in outlier detection will change.
#'
#' @export
#'
#' @return An exploratory data analysis PDF report.
#' @author Pol Castellano-Escuder
#'
#' @importFrom rmarkdown render
#' @import knitr
#' @import patchwork
#' @importFrom crayon red
#' @importFrom clisymbols symbol
#' 
#' @examples 
#' library(POMA)
#' data("st000284")
#' 
#' PomaEDA(st000284)
PomaEDA <- function(data, # nocov start
                    imputation = "knn",
                    normalization = "log_pareto",
                    clean_outliers = TRUE,
                    coeff_outliers = 1.5){
  
  if (missing(data)) {
    stop(crayon::red(clisymbols::symbol$cross, "data argument is empty!"))
  }
  if(!is(data)[1] == "MSnSet"){
    stop(paste0(crayon::red(clisymbols::symbol$cross, "data is not a MSnSet object."), 
                " \nSee POMA::PomaMSnSetClass or MSnbase::MSnSet"))
  }
  if (!(imputation %in% c("none", "half_min", "median", "mean", "min", "knn"))) {
    stop(crayon::red(clisymbols::symbol$cross, "Incorrect value for imputation argument!"))
  }
  if (!(normalization %in% c("none", "auto_scaling", "level_scaling", "log_scaling",
                      "log_transformation", "vast_scaling", "log_pareto"))) {
    stop(crayon::red(clisymbols::symbol$cross, "Incorrect value for normalization argument!"))
  }
  
  rmarkdown::render("R/POMA_EDA_report.Rmd", "pdf_document")
  
} # nocov end

