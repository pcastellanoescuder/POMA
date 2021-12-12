
#' Logistic Regression Model Odds Ratios
#'
#' @description PomaOddsRatio() calculates the Odds Ratios for each feature from a logistic regression model using the binary outcome (group/type must be a binary factor) as a dependent variable.
#'
#' @param data A SummarizedExperiment object. First `colData` column must be the subject group/type.
#' @param feature_name A vector with the name/s of feature/s that will be used to fit the model. If it's NULL (default), all variables will be included in the model.
#' @param covariates Logical that indicates if covariates will be included in logistic regression model. Default is `FALSE`.
#' @param covs Character vector indicating the name of `colData` columns that will be included as covariates. Default is NULL (all variables).
#' @param showCI Logical that indicates if the 95% confidence intervals will be plotted. Default is `TRUE`.
#'
#' @export
#'
#' @return A data frame with the Odds Ratios for all features with their 95% confidence intervals and a ggplot2 object.
#' @author Pol Castellano-Escuder
#'
#' @import ggplot2
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr mutate rename desc arrange mutate_at vars select_at matches as_tibble
#' @importFrom tidyr drop_na
#' @importFrom magrittr %>%
#' @importFrom SummarizedExperiment assay colData
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
  if(!is(data[1], "SummarizedExperiment")){
    stop("data is not a SummarizedExperiment object. \nSee POMA::PomaSummarizedExperiment or SummarizedExperiment::SummarizedExperiment")
  }
  if (!is.null(feature_name)) {
    if(!isTRUE(all(feature_name %in% rownames(SummarizedExperiment::assay(data))))){
      stop("At least one feature name not found...")
    }
  }
  if (length(levels(as.factor(SummarizedExperiment::colData(data)[,1]))) > 2) {
    stop("Your data have more than two groups!")
  }
  if(covariates & ncol(colData(data)) == 1){
    stop("Seems there aren't covariates in your data...")
  }

  e <- data.frame(t(SummarizedExperiment::assay(data)))
  target <- SummarizedExperiment::colData(data) %>% 
    as.data.frame() %>% 
    dplyr::rename(Group = 1)

  if(!is.null(feature_name)){
    e <- e %>% 
      dplyr::select_at(vars(matches(feature_name)))
    }

  if(covariates) {
    if(is.null(covs)){
      data_or <- cbind(target, e) %>% 
        as.data.frame() %>% 
        mutate(Group = as.factor(ifelse(Group == names(table(target$Group))[1], 0, 1))) %>% 
        mutate_at(vars(-Group), as.numeric)
    } 
    else{
      covariates <- target %>%
        as.data.frame() %>%
        dplyr::select_at(vars(matches("Group") | matches(covs)))
      
      data_or <- cbind(covariates, e) %>% 
        as.data.frame() %>% 
        mutate(Group = as.factor(ifelse(Group == names(table(target$Group))[1], 0, 1))) %>% 
        mutate_at(vars(-Group), as.numeric)
    }
  }
  else{
    data_or <- cbind(Group = target$Group, e) %>% 
      mutate(Group = as.factor(ifelse(Group == names(table(target$Group))[1], 0, 1))) %>% 
      mutate_at(vars(-Group), as.numeric)
  }
  
  suppressWarnings({
    logit_model <- glm(Group ~ 0 + ., family = binomial(link = 'logit'), data = data_or)
  })  

  odds <- data.frame(exp(cbind(OddsRatio = coef(logit_model), confint.default(logit_model, level = 0.95)))) %>%
    drop_na() %>%
    rownames_to_column("feature") %>%
    rename(CILow = `X2.5..`,
           CIHigh = `X97.5..`) %>%
    arrange(desc(OddsRatio)) %>% 
    dplyr::as_tibble()

  #### Plot

  ORPlot <- ggplot(odds, aes(x = OddsRatio, y = feature)) +
    geom_vline(aes(xintercept = 1), size = .25, linetype = "dashed") +
    {if(showCI)geom_errorbarh(aes(xmax = CIHigh, xmin = CILow), size = .5, height = .1, color = "gray50")} +
    geom_point(size = 3, color = "orange") +
    ylab("") +
    xlab("Odds ratio") +
    theme_bw()

  return(list(OddsRatioTable = odds, 
              OddsRatioPlot = ORPlot))

}

