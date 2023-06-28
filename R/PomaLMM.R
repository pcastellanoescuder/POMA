
#' Linear Mixed Models
#' 
#' @description `PomaLMM` performs linear mixed models on a `SummarizedExperiment` object.
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
#' # Perform linear mixed model with all features
#' st000284 %>% 
#' PomaLMM()
#' 
#' # Perform linear mixed model with two features
#' st000284 %>% 
#' PomaLMM(x = c("x1_methyladenosine", "x2_deoxyuridine"))
PomaLMM <- function(data,
                    x = NULL,
                    y = NULL,
                    adjust = "fdr"){
  
  if (!is(data, "SummarizedExperiment")){
    stop("data is not a SummarizedExperiment object. \nSee POMA::PomaCreateObject or SummarizedExperiment::SummarizedExperiment")
  }
  if (ncol(SummarizedExperiment::colData(data)) == 0) {
    stop("metadata file required")
  }
  if (is.null(y)) {
    y <- colnames(SummarizedExperiment::colData(data))
  }
  
  independent_variables <- SummarizedExperiment::colData(data) %>% 
    as.data.frame() %>% 
    dplyr::select(dplyr::all_of(y))
  
  fixed_effects <- independent_variables %>% 
    as.data.frame() %>% 
    dplyr::select_if(is.numeric)
  
  random_effects <- independent_variables %>% 
    as.data.frame() %>% 
    dplyr::select_if(is.factor)
  
  if (!(adjust %in% c("fdr", "holm", "hochberg", "hommel", "bonferroni", "BH", "BY"))) {
    stop("Incorrect value for adjust argument")
  }
  
  to_lmm <- data.frame(independent_variables, t(SummarizedExperiment::assay(data)))
  
  fixed_effects_names <- ifelse(ncol(fixed_effects) > 0,
                                paste0(colnames(fixed_effects), collapse = " + "),
                                NULL)
  
  random_effects_names <- ifelse(ncol(random_effects) > 0,
                                 paste0("(1|", colnames(random_effects), ")", collapse = " + "),
                                 NULL)
  
  formula_lmm <- paste0(fixed_effects_names, " + ", random_effects_names)
  formula_lmm <- as.formula(paste0(colnames(to_lmm)[1], " ~ ", formula_lmm))
  
  fit <- lme4::lmer(formula_lmm, independent_variables, REML = FALSE, data = to_lmm)
  
  # to_linear_model <- data.frame(dependent_variable, t(SummarizedExperiment::assay(data)))
  # 
  # if (is.null(x)) {
  #   res_lm <- broom::tidy(rlang::inject(lm(!!y_name ~ 0 + ., data = to_linear_model))) %>% 
  #     dplyr::mutate(adj_pvalue = p.adjust(p.value, method = adjust)) %>% 
  #     dplyr::select(feature = term, estimate, std_err = std.error, statistic, pvalue = p.value, adj_pvalue) %>% 
  #     dplyr::arrange(pvalue) %>% 
  #     dplyr::as_tibble()
  #   
  # } else {
  #   predictors <- to_linear_model %>%
  #     dplyr::select(dplyr::all_of(x))
  #   
  #   predictors <- data.frame(dependent_variable, predictors)
  #   
  #   res_lm <- broom::tidy(rlang::inject(lm(!!y_name ~ 0 + ., data = predictors))) %>% 
  #     dplyr::mutate(adj_pvalue = p.adjust(p.value, method = adjust)) %>% 
  #     dplyr::select(feature = term, estimate, std_err = std.error, statistic, pvalue = p.value, adj_pvalue) %>% 
  #     dplyr::arrange(pvalue) %>% 
  #     dplyr::as_tibble()
  # }
  # return(res_lm)
}

