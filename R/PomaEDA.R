
#' Automatic Exploratory Data Analysis PDF Report
#'
#' @description This function automatically generates a PDF report with different exploratory plots and tables from an SummarizedExperiment object.
#'
#' @param data A SummarizedExperiment object. First `colData` column must be the subject group/type.
#' @param imputation Imputation method. Options are "none", "half_min", "median", "mean", "min" and "knn" (default). If "none", all missing values will be replaced by zero.
#' @param normalization Normalization method. Options are "none", "auto_scaling", "level_scaling", "log_scaling", "log_transformation", "vast_scaling" and "log_pareto" (default).
#' @param clean_outliers Logical. If it's set to TRUE, outliers will be removed from EDA.
#' @param coeff_outliers This value corresponds to the classical 1.5 in \eqn{Q3 + 1.5*IQR} formula to detect outliers. By changing this value, the permissiveness in outlier detection will change.
#' @param username This name will be included as a report subtitle.
#' 
#' @export
#'
#' @return An exploratory data analysis PDF report.
#' @author Pol Castellano-Escuder
#'
#' @importFrom rmarkdown render
#' @import knitr
PomaEDA <- function(data, # nocov start
                    imputation = "knn",
                    normalization = "log_pareto",
                    clean_outliers = TRUE,
                    coeff_outliers = 1.5,
                    username = "Username"){
  
  if (missing(data)) {
    stop("data argument is empty!")
  }
  if(!is(data[1], "SummarizedExperiment")){
    stop("data is not a SummarizedExperiment object. \nSee POMA::PomaSummarizedExperiment or SummarizedExperiment::SummarizedExperiment")
  }
  if (!(imputation %in% c("none", "half_min", "median", "mean", "min", "knn"))) {
    stop("Incorrect value for imputation argument!")
  }
  if (!(normalization %in% c("none", "auto_scaling", "level_scaling", "log_scaling",
                      "log_transformation", "vast_scaling", "log_pareto"))) {
    stop("Incorrect value for normalization argument!")
  }
  
  rmarkdown::render(system.file("rmd", "POMA_EDA_report.Rmd", package = "POMA"), "pdf_document")
  
} # nocov end

