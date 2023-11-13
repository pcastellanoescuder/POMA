
#' Logistic Regression Model Odds Ratios
#'
#' @description `PomaOddsRatio` calculates the Odds Ratios for each feature from a logistic regression model using the binary outcome (group/type must be a binary factor) as a dependent variable.
#'
#' @param data A `SummarizedExperiment` object.
#' @param feature_name Character vector. Indicates the name/s of feature/s that will be used to fit the model. If it's NULL (default), all variables will be included in the model.
#' @param covs Character vector. Indicates the names of `colData` columns to be included as covariates. Default is NULL (no covariates).
#' @param show_ci Logical. Indicates if the 95% confidence intervals will be plotted. Default is `TRUE`.
#'
#' @export
#'
#' @return A `list` with results including plots and tables.
#' @author Pol Castellano-Escuder
#'
#' @importFrom magrittr %>%
#' 
#' @examples 
#' data("st000336")
#' 
#' st000336 %>% 
#'   PomaImpute() %>%
#'   PomaNorm() %>%
#'   PomaOddsRatio(feature_name = c("glutamic_acid", "glutamine", "glycine", "histidine"))
PomaOddsRatio <- function(data,
                          feature_name = NULL,
                          covs = NULL,
                          show_ci = TRUE) {

  if (!is(data, "SummarizedExperiment")){
    stop("data is not a SummarizedExperiment object. \nSee POMA::PomaCreateObject or SummarizedExperiment::SummarizedExperiment")
  }
  if (ncol(SummarizedExperiment::colData(data)) == 0) {
    stop("metadata file required")
  }
  if (!is.factor(SummarizedExperiment::colData(data)[,1])) {
    stop("Grouping factor must be a factor (first column of the metadata file)")
    if (length(table(SummarizedExperiment::colData(data)[,1])[table(SummarizedExperiment::colData(data)[,1]) != 0]) != 2) {
      stop("Grouping factor must have exactly 2 levels (first column of the metadata file)")
    }
  }
  if (!is.null(feature_name)) {
    if(!any(feature_name %in% rownames(SummarizedExperiment::assay(data)))) {
      stop("None of the specified features found")
    }
    if(!all(feature_name %in% rownames(SummarizedExperiment::assay(data)))){
      warning(paste0("Feature/s ",
                     paste0(feature_name[!feature_name %in% rownames(SummarizedExperiment::assay(data))], collapse = ", "),
                     " not found"))
    }
  }

  to_oddsratio <- data.frame(t(SummarizedExperiment::assay(data)))
  target <- SummarizedExperiment::colData(data) %>% 
    as.data.frame() %>% 
    dplyr::rename(group = 1)
  
  # feature names
  if (!is.null(feature_name)){
    to_oddsratio <- to_oddsratio[, colnames(to_oddsratio) %in% feature_name]
    
    if (is(to_oddsratio, "numeric")) {
      to_oddsratio <- data.frame(to_oddsratio)
      colnames(to_oddsratio) <- feature_name[feature_name %in% rownames(SummarizedExperiment::assay(data))]
    }
  }
  
  # covariates
  if (!is.null(covs)) {
    covariates <- target %>%
      dplyr::select_at(dplyr::vars(dplyr::matches("group") | dplyr::matches(covs)))
    
  } else {
    covariates <- target %>%
      dplyr::select(group)
  }
  
  data_or <- cbind(covariates, to_oddsratio) %>% 
    as.data.frame() %>% 
    dplyr::mutate(group = as.factor(ifelse(group == names(table(target$group))[1], 0, 1))) %>% 
    dplyr::mutate_at(dplyr::vars(-group), as.numeric)
  
  # Logistic model
  logit_model <- stats::glm(group ~ 0 + ., 
                            family = stats::binomial(link = 'logit'), 
                            data = data_or)

  odds <- data.frame(exp(cbind(OddsRatio = stats::coef(logit_model), 
                               stats::confint.default(logit_model, level = 0.95)))) %>%
    tidyr::drop_na() %>%
    tibble::rownames_to_column("feature") %>%
    dplyr::rename(lwr = 3, upr = 4) %>%
    dplyr::arrange(dplyr::desc(OddsRatio)) %>% 
    dplyr::as_tibble()

  ORPlot <- ggplot2::ggplot(odds, ggplot2::aes(x = OddsRatio, y = feature)) +
    ggplot2::geom_vline(xintercept = 1, size = .25, linetype = "dashed") +
    {if(show_ci) ggplot2::geom_errorbarh(ggplot2::aes(xmax = upr, xmin = lwr), size = .5, height = .1, color = "gray50")} +
    ggplot2::geom_point(size = 3, pch = 21, fill = "orange") +
    ggplot2::labs(x = "Odds Ratio",
                  y = NULL) +
    POMA::theme_poma()

  return(list(odds_ratio_table = odds, 
              odds_ratio_plot = ORPlot))
}

