
#' Normalization Plot
#'
#' @description PomaNormPlot() generates a boxplot of not normalized and normalized MS data. This plot can help in the comparison between pre and post normalized data and in the "validation" of the normalization process.
#'
#' @param data A MSnSet object. First `pData` column must be the suject group/type.
#' @param group Groupping factor for the plot. Options are c("subjects", "features"). If the user select "subject", the boxplot will be created for each subject. If the selection is "features", the boxplot will be created for each feature
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
PomaNormPlot <- function(data, group = c("subjects", "features")){

  if (!(group %in% c("subjects", "features"))) {
    stop(crayon::red(clisymbols::symbol$cross, "Incorrect value for group argument!"))
  }
  if (missing(group)) {
    group <- "subjects"
    warning("group argument is empty! subjects will be used")
  }

  e <- t(Biobase::exprs(data))
  pData <- Biobase::pData(data)

  data <- cbind(pData, e)
  data <- tibble::rownames_to_column(data, "ID")

  colnames(data)[2] <- c("Group")
  data <- data %>% mutate(ID = as.character(ID))

  if(group == "subjects"){

    data %>%
      reshape2::melt() %>%
      group_by(ID) %>%
      ggplot(aes(ID, value, color = Group)) +
      geom_boxplot() +
      geom_jitter() +
      theme_minimal() +
      ggtitle("Normalization Plot by Subjects") +
      xlab("") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))

  }
  else{

    data %>%
      dplyr::select(-ID) %>%
      reshape2::melt() %>%
      group_by(Group) %>%
      ggplot(aes(variable, value, color = Group)) +
      geom_boxplot() +
      geom_jitter() +
      theme_minimal() +
      xlab("") +
      ggtitle("Normalization Plot by Features") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))

  }

}

