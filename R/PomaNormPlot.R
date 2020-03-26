
#' Normalization Plot
#'
#' @description PomaNormPlot() generates a boxplot of not normalized and normalized MS data. This plot can help in the comparison between pre and post normalized data and in the "validation" of the normalization process.
#'
#' @param data A MSnSet object. First `pData` column must be the subject group/type.
#' @param group Groupping factor for the plot. Options are c("samples", "features"). Option "samples" (default) will create a boxplot for each sample and option "features" will create a boxplot of each variable.
#' @param jitter Logical. If it's TRUE (default), the boxplot will show all points.
#'
#' @export
#'
#' @return A ggplot2 object.
#' @author Pol Castellano-Escuder
#'
#' @import ggplot2
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr select group_by
#' @importFrom magrittr %>%
#' @importFrom reshape2 melt
#' @importFrom crayon red
#' @importFrom clisymbols symbol
#' @importFrom Biobase varLabels pData exprs
PomaNormPlot <- function(data, group = c("samples", "features"), jitter = TRUE){

  if (!(group %in% c("samples", "features"))) {
    stop(crayon::red(clisymbols::symbol$cross, "Incorrect value for group argument!"))
  }
  if (missing(group)) {
    group <- "samples"
    warning("group argument is empty! samples will be used")
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
      ggplot(aes(ID, value, color = Group)) +
      geom_boxplot() +
      {if(jitter)geom_jitter()} +
      theme_bw() +
      xlab("Samples") +
      ylab("Value") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    }

  else {

    data %>%
      dplyr::select(-ID) %>%
      reshape2::melt() %>%
      group_by(Group) %>%
      ggplot(aes(variable, value, color = Group)) +
      geom_boxplot() +
      {if(jitter)geom_jitter()} +
      theme_bw() +
      xlab("Features") +
      ylab("Value") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    }
  }

