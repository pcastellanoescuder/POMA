
#' Principal Components Regression
#'
#' @description Fits a Linear Model using the indicated variable as the dependent variable (outcome) and the indicated number of Principal Components as independent variables.
#'
#' @param data A SummarizedExperiment object.
#' @param n_components The number of Principal Components used in the PCR model. Defaults is 5.
#' @param scale A logical value indicating whether the variables should be scaled to have unit variance before the analysis takes place. 
#' @param center A logical value indicating whether the variables should be shifted to be zero centered. 
#' @param outcome Character string indicating the variable name in `colData` to be used as the outcome.
#' @param intercept A logical value indicating whether intercept should be included in the LM. Default is TRUE.
#' 
#' @export
#'
#' @return A list with PCR results.
#' @author Pol Castellano-Escuder
#'
#' @importFrom magrittr %>% %<>%
#' 
#' @examples 
#' data("st000284")
#' 
#' st000284 %>%
#'   PomaPCR(outcome = "age_at_consent")
#' 
#' data("st000336")
#'   
#' st000336 %>% 
#' PomaImpute() %>% 
#' PomaPCR(n_components = 2, 
#'         outcome = "steroids", 
#'         intercept = FALSE)
PomaPCR <- function(data,
                    n_components = 5,
                    scale = TRUE,
                    center = TRUE,
                    outcome = NULL,
                    intercept = TRUE) {
  
  if (missing(data)) {
    stop("data argument is empty!")
  }
  if(!is(data, "SummarizedExperiment")){
    stop("data is not a SummarizedExperiment object. \nSee POMA::PomaSummarizedExperiment or SummarizedExperiment::SummarizedExperiment")
  }
  if(is.null(outcome)){
    stop("Select an outcome variable")
  }
  if(!outcome %in% colnames(SummarizedExperiment::colData(data))) {
    stop(paste0("The variable ", outcome, " is not in your target data."))
  }
  
  assay <- t(SummarizedExperiment::assay(data))
  target <- SummarizedExperiment::colData(data)
  
  pca_res <- prcomp(assay, scale. = scale, center = center)
  
  pca_loadings <- pca_res$rotation[, 1:n_components] %>% 
    dplyr::as_tibble() %>% 
    dplyr::mutate(feature = rownames(data)) %>% 
    dplyr::relocate(feature, dplyr::everything()) %>% 
    dplyr::arrange(dplyr::desc(abs(PC1)))
  
  pca_res <- pca_res$x[, 1:n_components]

  outcome_df <- target %>%
    as.data.frame() %>%
    dplyr::select_at(dplyr::vars(outcome = dplyr::matches(outcome)))
  
  data_pcr <- cbind(outcome_df, pca_res)
  
  if(is(outcome_df$outcome, "integer") | is(outcome_df$outcome, "numeric")) {
    data_pcr %<>%
      dplyr::mutate(outcome = as.numeric(outcome))
    
    if(intercept) {
      res_pcr <- lm(outcome ~ ., 
                    data = data_pcr)  
    } else {
      res_pcr <- lm(outcome ~ . - 1, 
                    data = data_pcr) 
    }
  }
  
  else {
    stop("The outcome variable is expected to be numeric")
  }

  return(list(summary = broom::glance(res_pcr),
              coefficients = broom::tidy(res_pcr),
              loadings = pca_loadings)
         )
  
}

