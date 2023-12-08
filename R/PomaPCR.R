
#' Principal Components Regression
#'
#' @description `PomaPCR` performs Principal Components Regression.
#'
#' @param data A `SummarizedExperiment` object.
#' @param center Logical. Indicates whether the variables should be shifted to be zero centered. Default is TRUE.
#' @param scale Logical. Indicates whether the variables should be scaled to have unit variance before the analysis takes place. Default is TRUE.
#' @param ncomp Numeric. Indicates the number of principal components used as predictors in the model. Default is 2.
#' @param y Character. Indicates the name of `colData` columns to be used as dependent variable. If it's set to NULL, the first numeric variable in `colData` will be used as the dependent variable.
#' @param adjust Character. Multiple comparisons correction method to adjust p-values. Available options are: "fdr" (false discovery rate), "holm", "hochberg", "hommel", "bonferroni", "BH" (Benjamini-Hochberg), and "BY" (Benjamini-Yekutieli).
#' 
#' @export
#'
#' @return A `tibble` with the results.
#' @author Pol Castellano-Escuder
#'
#' @importFrom magrittr %>% %<>%
#' 
#' @examples 
#' data("st000284")
#' 
#' # PCR with 2 components
#' st000284 %>%
#'   PomaPCR(y = "age_at_consent")
#'   
#' # PCR with 20 components
#' st000284 %>%
#'   PomaPCR(ncomp = 20)
PomaPCR <- function(data,
                    center = TRUE,
                    scale = TRUE,
                    ncomp = 2,
                    y = NULL,
                    adjust = "fdr") {
  
  if(!is(data, "SummarizedExperiment")){
    stop("data is not a SummarizedExperiment object. \nSee POMA::PomaCreateObject or SummarizedExperiment::SummarizedExperiment")
  }
  if (ncol(SummarizedExperiment::colData(data)) == 0) {
    stop("metadata file required")
  }
  
  dependent_variable <- SummarizedExperiment::colData(data) %>% 
    as.data.frame() %>% 
    dplyr::select_if(is.numeric)
  
  if (ncol(dependent_variable) == 0) {
    stop("No numeric variables to be used as dependent variable in metadata file")
  }
  if (is.null(y)) {
    y <- colnames(dependent_variable)[1]
  }
  
  y_name <- rlang::sym(y[1])
  
  dependent_variable <- dependent_variable %>% 
    dplyr::select(dplyr::all_of(y_name))
  
  if (!(adjust %in% c("fdr", "holm", "hochberg", "hommel", "bonferroni", "BH", "BY"))) {
    stop("Incorrect value for adjust argument")
  }

  pca_res <- POMA::PomaPCA(data, center = center, scale = scale, ncomp = ncomp)$factors %>% 
    dplyr::select(PC1:paste0("PC", ncomp))
  
  to_pcr <- data.frame(dependent_variable, pca_res)

  res_pcr <- broom::tidy(rlang::inject(lm(!!y_name ~ 0 + ., data = to_pcr))) %>% 
    dplyr::mutate(adj_pvalue = p.adjust(p.value, method = adjust)) %>% 
    dplyr::select(component = term, estimate, std_err = std.error, statistic, pvalue = p.value, adj_pvalue) %>% 
    dplyr::arrange(pvalue) %>% 
    dplyr::as_tibble()

  return(res_pcr)
}

