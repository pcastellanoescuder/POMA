
#' Batch Correction
#'
#' @description `PomaBatch` performs batch correction on a `SummarizedExperiment` object given a batch factor variable.
#'
#' @param data A `SummarizedExperiment` object.
#' @param batch Character. The name of the column in `colData` that contains the batch information.
#' @param mod Character vector. Indicates the names of `colData` columns to be included as covariates. Default is NULL (no covariates).
#'
#' @export
#'
#' @return A `SummarizedExperiment` object with batch-corrected data.
#' @references Leek JT, Johnson WE, Parker HS, Fertig EJ, Jaffe AE, Zhang Y, Storey JD, Torres LC (2023). sva: Surrogate Variable Analysis. doi:10.18129/B9.bioc.sva <https://doi.org/10.18129/B9.bioc.sva>
#' @author Pol Castellano-Escuder
#'
#' @importFrom magrittr %>%
#' 
#' @examples 
#' # Output is a batch corrected SummarizedExperiment object
#' data <- POMA::st000284 # Example SummarizedExperiment object included in POMA
#' 
#' data %>%
#'   PomaBatch(batch = "gender")
PomaBatch <- function(data, 
                      batch, 
                      mod = NULL) {
  
  if (!is(data, "SummarizedExperiment")) {
    stop("data is not a SummarizedExperiment object. \nSee POMA::PomaCreateObject or SummarizedExperiment::SummarizedExperiment")
  }
  if (!batch %in% names(SummarizedExperiment::colData(data))) {
    stop("Specified batch column not found in colData")
  }
  
  batch_info <- as.factor(SummarizedExperiment::colData(data)[, batch])
  to_batch <- SummarizedExperiment::assay(data)
  
  if (!is.null(mod)) {
    covariates <- SummarizedExperiment::colData(data) %>%
      as.data.frame() %>%
      dplyr::select_at(dplyr::vars(dplyr::matches(mod)))
    
    mod_matrix <- stats::model.matrix(~ 0 + ., data = covariates)
  } else {
    mod_matrix <- NULL
  }

  corrected_data <- sva::ComBat(dat = to_batch, batch = batch_info, mod = mod_matrix)
  
  data <- SummarizedExperiment::SummarizedExperiment(assays = corrected_data, 
                                                     colData = SummarizedExperiment::colData(data))
  
  if (validObject(data))
    return(data)
}

