
#' Linear Mixed Models
#' 
#' @description `PomaLMM` performs linear mixed models on a `SummarizedExperiment` object.
#'
#' @param data A `SummarizedExperiment` object.
#' @param x Character vector. Indicates the names of `colData` columns to be used as random and fixed effects (independent variables). If it's set to NULL (default), all variables in `colData` will be used.
#' @param y Character vector. Indicates the names of dependent variables. If it's NULL (default), all features will be used.
#' @param adjust Character. Multiple comparisons correction method to adjust p-values. Available options are: "fdr" (false discovery rate), "holm", "hochberg", "hommel", "bonferroni", "BH" (Benjamini-Hochberg), and "BY" (Benjamini-Yekutieli).
#' @param clean_plot Logical. Indicates if remove intercept and linear mixed model residues boxplots from the plot. Defasult is FALSE. 
#'
#' @export
#'
#' @return A `list` with results including plots and tables. Table values indicate the percentage variance explained per variable.
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
#' PomaLMM(y = c("x1_methyladenosine", "x1_methylhistamine"))
#' 
#' # Perform linear mixed model with one random effect
#' st000284 %>% 
#' PomaLMM(x = "smoking_condition")
#' 
#' # Perform linear mixed model with two random effects and two features
#' st000284 %>% 
#' PomaLMM(x = c("smoking_condition", "gender"),
#'         y = c("x1_methyladenosine", "x1_methylhistamine"))
#'         
#' # Perform linear mixed model with no random effects and two features, therefore, a linear model will be fitted
#' st000284 %>% 
#' PomaLMM(x = "age_at_consent", # Numerical, i.e., fixed effect
#'         y = c("x1_methyladenosine", "x1_methylhistamine"))
#'         
#' # Perform linear mixed model with no random effects and all features, therefore, a linear model will be fitted
#' st000284 %>% 
#' PomaLMM(x = "age_at_consent") # Numerical i.e., fixed effect
PomaLMM <- function(data,
                    x = NULL,
                    y = NULL,
                    adjust = "fdr",
                    clean_plot = FALSE) {
  
  if (!is(data, "SummarizedExperiment")){
    stop("data is not a SummarizedExperiment object. \nSee POMA::PomaCreateObject or SummarizedExperiment::SummarizedExperiment")
  }
  if (ncol(SummarizedExperiment::colData(data)) == 0) {
    stop("metadata file required")
  }
  if (is.null(x)) {
    x <- colnames(SummarizedExperiment::colData(data))
  }
  
  independent_variables <- SummarizedExperiment::colData(data) %>% 
    as.data.frame() %>% 
    dplyr::select(dplyr::all_of(x))
  
  fixed_effects <- independent_variables %>% 
    as.data.frame() %>% 
    dplyr::select_if(is.numeric)
  
  random_effects <- independent_variables %>% 
    as.data.frame() %>% 
    dplyr::select_if(is.factor)
  
  if (ncol(random_effects) == 0) {
    message("No random effects. A linear model will be fitted.")
    
    res_lm <- PomaLM(data = data, 
                     x = y,
                     y = colnames(fixed_effects)[1],
                     adjust = adjust)
    
    return(res_lm)
  }
  
  if (!(adjust %in% c("fdr", "holm", "hochberg", "hommel", "bonferroni", "BH", "BY"))) {
    stop("Incorrect value for adjust argument")
  }
  
  to_lmm <- t(SummarizedExperiment::assay(data))

  if (!is.null(y)) {
    to_lmm <- to_lmm[, colnames(to_lmm) %in% y]
  }
  
  formula_str <- paste0("response ~ ", paste(names(fixed_effects), collapse = " + "))
  if (ncol(random_effects) > 0) {
    formula_str <- paste0(formula_str, " + (1 | ", paste(names(random_effects), collapse = ") + (1 | "), ")")
  }
  
  lmm_fun <- function(x) {
    model <- lme4::lmer(formula_str, data = data.frame(response = x, fixed_effects, random_effects), REML = FALSE)
    
    random_effects_variances <- data.frame(lme4::VarCorr(model))[, c(1,4)] %>% 
      dplyr::as_tibble()
    
    fixed_effects_variances <- data.frame(vcov = diag(as.matrix(vcov(model)))) %>% 
      tibble::rownames_to_column("grp") %>% 
      dplyr::as_tibble()
    
    total_vars <- fixed_effects_variances %>% 
      dplyr::bind_rows(random_effects_variances) %>% 
      dplyr::mutate(variance_percent = 100*(vcov/sum(vcov))) %>% 
      dplyr::select(-vcov) %>% 
      tidyr::pivot_wider(names_from = grp, values_from = variance_percent)
    
    return(total_vars)
  }
  
  # variances
  suppressMessages({
    res_lmm <- dplyr::bind_rows(apply(to_lmm, 2, function(x){lmm_fun(x)})) %>% 
      dplyr::mutate(feature = colnames(to_lmm)) %>% 
      dplyr::select(feature, dplyr::everything()) %>% 
      dplyr::as_tibble()
  })
  
  # variances plot
  plot_data <- res_lmm %>% 
    tidyr::pivot_longer(cols = -feature)
  
  if (clean_plot) {
    plot_data <- plot_data %>% 
      dplyr::filter(!name %in% c("(Intercept)", "Residual"))
  }
  
  variances_plot <- plot_data %>%
    ggplot2::ggplot(ggplot2::aes(name, value)) +
    ggplot2::geom_boxplot(ggplot2::aes(fill = name), show.legend = FALSE) +
    ggplot2::labs(x = NULL,
                  y = "Variance Explained (%)",
                  fill = NULL) +
    theme_poma(axis_x_rotate = TRUE) +
    scale_fill_poma_d()
    
  return(list(variances = res_lmm,
              variances_plot = variances_plot)
         )
}

