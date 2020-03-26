
#' Logistic Regression Model Odds Ratios
#'
#' @description PomaOddsRatio() calculates the Odds Ratios for each feature from a logistic regression model using the binary outcone (group/type must be a binary factor) as a dependent variable.
#'
#' @param data A MSnSet object. First `pData` column must be the subject group/type.
#' @param feature_name A vector with the name/s of feature/s to plot. If it's NULL (default), all variables will be plotted.
#'
#' @export
#'
#' @return A data frame with the Odds Ratios for all features with their 95% confidence intervals and a ggplot2 object.
#' @author Pol Castellano-Escuder
#'
#' @import ggplot2
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr mutate rename desc arrange
#' @importFrom tidyr drop_na
#' @importFrom magrittr %>%
#' @importFrom crayon red
#' @importFrom clisymbols symbol
#' @importFrom Biobase pData exprs featureNames
PomaOddsRatio <- function(data, feature_name = NULL){

  if (!is.null(feature_name)) {
    if(!(feature_name %in% Biobase::featureNames(data))){
      stop(crayon::red(clisymbols::symbol$cross, "Feature name not found!"))
    }
  }

  e <- t(Biobase::exprs(data))
  pData <- Biobase::pData(data)[1]
  colnames(pData) <- "Group"

  data <- as.data.frame(cbind(Group = pData$Group, e))
  data <- data %>% mutate(Group = as.factor(ifelse(Group == levels(data$Group)[1], 0, 1)))

  data[,2:ncol(data)] <- as.data.frame(apply(data[,2:ncol(data)], 2, function(x){as.numeric(as.character(x))}))

  logit_model <- glm(Group ~ 0 +., family = binomial(link = 'logit'), data = data)

  odds <- data.frame(exp(cbind(OddsRatio = coef(logit_model), confint.default(logit_model, level = 0.95)))) %>%
    drop_na() %>%
    rownames_to_column("feature") %>%
    rename(CILow = `X2.5..`,
           CIHigh = `X97.5..`) %>%
    arrange(desc(OddsRatio))

  #### Plot

  if(!is.null(feature_name)){
    odds_mod <- odds %>% filter(feature %in% feature_name)
  }
  else{
    odds_mod <- odds
  }

  ORPlot <- ggplot(odds_mod, aes(x = OddsRatio, y = feature)) +
    geom_vline(aes(xintercept = 1), size = .25, linetype = "dashed") +
    geom_errorbarh(aes(xmax = CIHigh, xmin = CILow), size = .5, height = .2, color = "gray50") +
    geom_point(size = 3, color = "orange") +
    ylab("") +
    xlab("Odds ratio") +
    theme_bw()

  return(list(OddsRatioTable = odds, OddsRatioPlot = ORPlot))

}

