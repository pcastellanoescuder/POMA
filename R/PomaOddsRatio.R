
#' Logistic Regression Model Odds Ratios
#'
#' @description PomaOddsRatio() calculates the Odds Ratios for each feature from a logistic regression model using the binary outcome (group/type must be a binary factor) as a dependent variable.
#'
#' @param data A MSnSet object. First `pData` column must be the subject group/type.
#' @param feature_name A vector with the name/s of feature/s that will be used to fit the model. If it's NULL (default), all variables will be included in the model.
#' @param covariates Logical that indicates if covariates will be included in logistic regression model. Default is `FALSE`.
#' @param showCI Logical that indicates if the 95% confidence intervals will be plotted. Default is `TRUE`.
#'
#' @export
#'
#' @return A data frame with the Odds Ratios for all features with their 95% confidence intervals and a ggplot2 object.
#' @author Pol Castellano-Escuder
#'
#' @import ggplot2
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr mutate rename desc arrange bind_cols mutate_at vars
#' @importFrom tidyr drop_na
#' @importFrom magrittr %>%
#' @importFrom crayon red
#' @importFrom clisymbols symbol
#' @importFrom Biobase pData exprs featureNames
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
                          showCI = TRUE){

  if (missing(data)) {
    stop(crayon::red(clisymbols::symbol$cross, "data argument is empty!"))
  }
  if(!is(data[1], "MSnSet")){
    stop(paste0(crayon::red(clisymbols::symbol$cross, "data is not a MSnSet object."), 
                " \nSee POMA::PomaMSnSetClass or MSnbase::MSnSet"))
  }
  if (!is.null(feature_name)) {
    if(!isTRUE(all(feature_name %in% Biobase::featureNames(data)))){
      stop(crayon::red(clisymbols::symbol$cross, "At least one feature name not found..."))
    }
  }
  if(isTRUE(covariates) & ncol(pData(data)) == 1){
    stop(crayon::red(clisymbols::symbol$cross, "Seems that your data don't have covariates..."))
  }

  e <- data.frame(t(Biobase::exprs(data)))
  pData <- Biobase::pData(data)
  colnames(pData)[1] <- "Group"

  if(!is.null(feature_name)){
    e <- as.data.frame(e[, colnames(e) %in% feature_name])

    if(length(feature_name == 1)){
      colnames(e) <- feature_name 
    }

  } else {
    e <- as.data.frame(e)
  }

  if(covariates){
    data <- bind_cols(pData, e)
    data$Group <- as.factor(data$Group)

  } else {
    data <- bind_cols(Group = as.factor(pData$Group), e)
  }

  data <- data %>% 
    mutate(Group = as.factor(ifelse(Group == levels(data$Group)[1], 0, 1))) %>% 
    mutate_at(vars(-Group), as.numeric)
  
  suppressWarnings({
    logit_model <- glm(Group ~ 0 + ., family = binomial(link = 'logit'), data = data)
  })  

  odds <- data.frame(exp(cbind(OddsRatio = coef(logit_model), confint.default(logit_model, level = 0.95)))) %>%
    drop_na() %>%
    rownames_to_column("feature") %>%
    rename(CILow = `X2.5..`,
           CIHigh = `X97.5..`) %>%
    arrange(desc(OddsRatio))

  #### Plot

  ORPlot <- ggplot(odds, aes(x = OddsRatio, y = feature)) +
    geom_vline(aes(xintercept = 1), size = .25, linetype = "dashed") +
    {if(isTRUE(showCI))geom_errorbarh(aes(xmax = CIHigh, xmin = CILow), size = .5, height = .1, color = "gray50")} +
    geom_point(size = 3, color = "orange") +
    ylab("") +
    xlab("Odds ratio") +
    theme_bw()

  return(list(OddsRatioTable = odds, OddsRatioPlot = ORPlot))

}

