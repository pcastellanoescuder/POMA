
#' Automatic Exploratory Data Analysis HTML Report
#'
#' @description This function automatically generates a HTML report with different exploratory plots and tables from an `SummarizedExperiment` object.
#'
#' @param data A SummarizedExperiment object.
#' @param imputation Imputation method. See `?POMA::PomaImpute()`.
#' @param normalization Normalization method. See `?POMA::PomaNorm()`.
#' @param clean_outliers Logical. If it's set to TRUE, outliers will be removed from EDA.
#' @param coeff_outliers This value corresponds to the classical 1.5 in \eqn{Q3 + 1.5*IQR} formula to detect outliers. See `?POMA::PomaOutliers()`.
#' @param username Author name in the report.
#' @param institution Institution name in the report.
#' 
#' @export
#'
#' @return An exploratory data analysis HTML report.
#' @author Pol Castellano-Escuder
PomaEDA <- function(data, # nocov start
                    imputation = "knn",
                    normalization = "log_pareto",
                    clean_outliers = TRUE,
                    coeff_outliers = 1.5,
                    username = NULL,
                    institution = NULL){
  
  if (missing(data)) {
    stop("data argument is empty!")
  }
  if(!is(data, "SummarizedExperiment")){
    stop("data is not a SummarizedExperiment object. \nSee POMA::PomaSummarizedExperiment or SummarizedExperiment::SummarizedExperiment")
  }
  if (!(imputation %in% c("none", "half_min", "median", "mean", "min", "knn"))) {
    stop("Incorrect value for imputation argument!")
  }
  if (!(normalization %in% c("none", "auto_scaling", "level_scaling", "log_scaling",
                      "log_transformation", "vast_scaling", "log_pareto"))) {
    stop("Incorrect value for normalization argument!")
  }
  if(!is(SummarizedExperiment::colData(data)[,1], "character") &
     !is(SummarizedExperiment::colData(data)[,1], "factor")){
    stop("PomaEDA expects the first column of your target to be a factor or character. More report options coming soon...")
  }
  
  rmarkdown::render(system.file("rmd", "POMA_EDA_report.Rmd", package = "POMA"), "html_document")
  
} # nocov end

