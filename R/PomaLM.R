
#' Linear Models
#' 
#' @description `PomaLM` performs a linear model on a `SummarizedExperiment` object.
#'
#' @param data A `SummarizedExperiment` object.
#' @param x Character vector. Indicates the names of independent variables. If it's NULL (default), all features will be used.
#' @param y Character. Indicates the name of `colData` numeric columns to be used as dependent variable. If it's set to NULL, the first numeric variable in `colData` will be used as the dependent variable.
#' @param adjust Character. Multiple comparisons correction method to adjust p-values. Available options are: "fdr" (false discovery rate), "holm", "hochberg", "hommel", "bonferroni", "BH" (Benjamini-Hochberg), and "BY" (Benjamini-Yekutieli).
#'
#' @export
#'
#' @return A `list` with results including plots and tables. 
#' @author Pol Castellano-Escuder
#'
#' @importFrom magrittr %>%
#' 
#' @examples 
#' data <- POMA::st000284 %>% # Example SummarizedExperiment object included in POMA
#'   PomaImpute() %>% 
#'   PomaNorm()
#' 
#' ## Output is a list with objects `lm_table` (tibble) and `regression_plot` (ggplot2 object)
#' # Perform linear model with all features
#' data %>% 
#'   PomaLM()
#' 
#' # Perform linear model with two features
#' data %>% 
#'   PomaLM(x = c("x1_methyladenosine", "x2_deoxyuridine"))
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
  
  # regression plot
  regression_plot <- res_lm %>%
    ggplot2::ggplot(ggplot2::aes(x = estimate, y = reorder(feature, estimate), fill = pvalue)) +
    ggplot2::geom_vline(xintercept = 0, size = 0.25, linetype = "dashed") +
    ggplot2::geom_errorbarh(ggplot2::aes(xmax = estimate + std_err, xmin = estimate - std_err), 
                            size = 0.5, height = 0.1, color = "black") +
    ggplot2::geom_point(size = 3, pch = 21, alpha = 0.8) + 
    ggplot2::labs(x = "Coefficient",
                  y = NULL,
                  fill = "p-value") +
    theme_poma(legend_position = "right") +
    scale_fill_poma_c()
  
  return(list(lm_table = res_lm,
              regression_plot = regression_plot)
         )
}

