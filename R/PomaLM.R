
#' Linear Models
#' 
#' @description `PomaLM` performs a linear model on a `SummarizedExperiment` object.
#'
#' @param data A `SummarizedExperiment` object.
#' @param x Character vector. Indicates the names of independent variables. If it's null, all features will be used.
#' @param y Character. Indicates the name of `colData` columns to be used as dependent variable. If it's set to NULL, the first numeric variable in `colData` will be used as the dependent variable.
#' @param adjust Character. Multiple comparisons correction method to adjust p-values. Available options are: "fdr" (false discovery rate), "holm", "hochberg", "hommel", "bonferroni", "BH" (Benjamini-Hochberg), and "BY" (Benjamini-Yekutieli).
#'
#' @export
#'
#' @return A `tibble` with the results.
#' @author Pol Castellano-Escuder
#'
#' @importFrom magrittr %>%
#' 
#' @examples 
#' data("st000284")
#' 
#' # Perform linear model with all features
#' st000284 %>% 
#' PomaLM()
#' 
#' # Perform linear model with two features
#' st000284 %>% 
#' PomaLM(x = c("x1_methyladenosine", "x2_deoxyuridine"))
PomaLM <- function(data,
                   x = NULL,
                   y = NULL,
                   adjust = "fdr"){
  
  if (!is(data, "SummarizedExperiment")){
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
  
  to_linear_model <- data.frame(dependent_variable, t(SummarizedExperiment::assay(data)))
    
  if (is.null(x)) {
    res_lm <- broom::tidy(rlang::inject(lm(!!y_name ~ 0 + ., data = to_linear_model))) %>% 
      dplyr::mutate(adj_pvalue = p.adjust(p.value, method = adjust)) %>% 
      dplyr::select(feature = term, estimate, std_err = std.error, statistic, pvalue = p.value, adj_pvalue) %>% 
      dplyr::arrange(pvalue) %>% 
      dplyr::as_tibble()
    
  } else {
    predictors <- to_linear_model %>%
      dplyr::select(dplyr::all_of(x))
    
    predictors <- data.frame(dependent_variable, predictors)
   
    res_lm <- broom::tidy(rlang::inject(lm(!!y_name ~ 0 + ., data = predictors))) %>% 
      dplyr::mutate(adj_pvalue = p.adjust(p.value, method = adjust)) %>% 
      dplyr::select(feature = term, estimate, std_err = std.error, statistic, pvalue = p.value, adj_pvalue) %>% 
      dplyr::arrange(pvalue) %>% 
      dplyr::as_tibble()
  }
  return(res_lm)
}

