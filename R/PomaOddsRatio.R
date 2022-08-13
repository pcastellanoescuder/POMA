
#' Logistic Regression Model Odds Ratios
#'
#' @description PomaOddsRatio() calculates the Odds Ratios for each feature from a logistic regression model using the binary outcome (group/type must be a binary factor) as a dependent variable.
#'
#' @param data A SummarizedExperiment object.
#' @param feature_name A vector with the name/s of feature/s that will be used to fit the model. If it's NULL (default), all variables will be included in the model.
#' @param covariates Logical that indicates if covariates will be included in logistic regression model. Default is `FALSE`.
#' @param covs Character vector indicating the name of `colData` columns that will be included as covariates. Default is NULL (all variables).
#' @param showCI Logical that indicates if the 95% confidence intervals will be plotted. Default is `TRUE`.
#'
#' @export
#'
#' @return A tibble with the Odds Ratios for all features with their 95% confidence intervals and a ggplot2 object.
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
#'   PomaOddsRatio(feature_name = c("glutamic_acid", "glutamine", 
#'                                  "glycine", "histidine"))
PomaOddsRatio <- function(data,
                          feature_name = NULL,
                          covariates = FALSE,
                          covs = NULL,
                          showCI = TRUE){

  if (missing(data)) {
    stop("data argument is empty!")
  }
  if(!is(data, "SummarizedExperiment")){
    stop("data is not a SummarizedExperiment object. \nSee POMA::PomaSummarizedExperiment or SummarizedExperiment::SummarizedExperiment")
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
  if (length(levels(as.factor(SummarizedExperiment::colData(data)[,1]))) > 2) {
    stop("Your data have more than two groups!")
  }
  if(covariates & ncol(SummarizedExperiment::colData(data)) == 1){
    stop("No covariates found in colData")
  }

  e <- data.frame(t(SummarizedExperiment::assay(data)))
  target <- SummarizedExperiment::colData(data) %>% 
    as.data.frame() %>% 
    dplyr::rename(group = 1)

  if(!is.null(feature_name)){
    e <- e[, colnames(e) %in% feature_name]
    
    if(is(e, "numeric")) {
      e <- data.frame(e)
      colnames(e) <- feature_name[feature_name %in% rownames(SummarizedExperiment::assay(data))]
    }
  }

  if(covariates) {
    if(is.null(covs)){
      data_or <- cbind(target, e) %>% 
        as.data.frame() %>% 
        dplyr::mutate(group = as.factor(ifelse(group == names(table(target$group))[1], 0, 1))) %>% 
        dplyr::mutate_at(dplyr::vars(-group), as.numeric)
    } 
    else{
      covariates <- target %>%
        as.data.frame() %>%
        dplyr::select_at(dplyr::vars(matches("group") | dplyr::matches(covs)))
      
      data_or <- cbind(covariates, e) %>% 
        as.data.frame() %>% 
        dplyr::mutate(group = as.factor(ifelse(group == names(table(target$group))[1], 0, 1))) %>% 
        dplyr::mutate_at(dplyr::vars(-group), as.numeric)
    }
  }
  else{
    data_or <- cbind(group = target$group, e) %>% 
      dplyr::mutate(group = as.factor(ifelse(group == names(table(target$group))[1], 0, 1))) %>% 
      dplyr::mutate_at(dplyr::vars(-group), as.numeric)
  }
  
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
    {if(showCI)ggplot2::geom_errorbarh(ggplot2::aes(xmax = upr, xmin = lwr), size = .5, height = .1, color = "gray50")} +
    ggplot2::geom_point(size = 3, color = "orange") +
    ggplot2::labs(x = "Odds Ratio",
                  y = NULL) +
    ggplot2::theme_bw()

  return(list(OddsRatioTable = odds, 
              OddsRatioPlot = ORPlot))

}

