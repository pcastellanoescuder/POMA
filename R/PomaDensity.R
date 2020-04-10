
#' Distribution Plot
#'
#' @description PomaDensity() generates a density plot of not normalized and normalized MS data. This plot can help in the comparison between pre and post normalized data and in the "validation" of the normalization process.
#'
#' @param data A MSnSet object. First `pData` column must be the subject group/type.
#' @param group Groupping factor for the plot. Options are "samples" and "features". Option "samples" (default) will create a density plot for each group and option "features" will create a density plot of each variable.
#' @param feature_name A vector with the name/s of feature/s to plot. If it's NULL (default) a density plot of all variables will be created.
#'
#' @export
#'
#' @return A ggplot2 object.
#' @author Pol Castellano-Escuder
#'
#' @import ggplot2
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr select group_by filter
#' @importFrom magrittr %>%
#' @importFrom reshape2 melt
#' @importFrom crayon red
#' @importFrom clisymbols symbol
#' @importFrom Biobase varLabels pData exprs featureNames
PomaDensity <- function(data,
                        group = "samples",
                        feature_name = NULL){

  if (!(group %in% c("samples", "features"))) {
    stop(crayon::red(clisymbols::symbol$cross, "Incorrect value for group argument!"))
  }
  if (missing(group)) {
    warning("group argument is empty! samples will be used")
  }
  if (!is.null(feature_name)) {
    if(!(feature_name %in% Biobase::featureNames(data))){
      stop(crayon::red(clisymbols::symbol$cross, "Feature name not found!"))
    }
  }

  e <- t(Biobase::exprs(data))
  pData <- Biobase::pData(data)

  data <- cbind(pData, e)
  data <- tibble::rownames_to_column(data, "ID")

  colnames(data)[2] <- c("Group")
  data <- data %>% mutate(ID = as.character(ID))

  n_rem <- ncol(pData) - 1
  data <- data[, c(1:2, (3 + n_rem):ncol(data))]

  if(group == "samples"){

    data %>%
      reshape2::melt() %>%
      group_by(ID) %>%
      ggplot(aes(value, fill = Group)) +
      geom_density(alpha = 0.4) +
      xlab("Value") +
      ylab("Density") +
      theme_bw()

  } else {

    if (is.null(feature_name)){

      data %>%
        dplyr::select(-ID) %>%
        reshape2::melt() %>%
        group_by(Group) %>%
        ggplot(aes(value, fill = variable)) +
        geom_density(alpha = 0.4) +
        theme_bw() +
        xlab("Value") +
        ylab("Density") +
        theme(legend.position = "none")

    } else {

      data %>%
        dplyr::select(-ID) %>%
        reshape2::melt() %>%
        group_by(Group) %>%
        filter(variable %in% feature_name) %>%
        ggplot(aes(value, fill = variable)) +
        geom_density(alpha = 0.4) +
        theme_bw() +
        xlab("Value") +
        ylab("Density")

    }
  }
}

