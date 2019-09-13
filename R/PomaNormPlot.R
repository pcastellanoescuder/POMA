
#' Normalization Plot
#'
#' @description PomaNormPlot() generates a boxplot of not normalized and normalized metabolomic data. This plot can help in the comparison between pre and post normalized data and in the "validation" of the normalization process.
#'
#' @param data A data frame with metabolites. First column must be the subject ID and second column must be a factor with the subject group.
#' @param group Groupping factor for the plot. Options are c("subjects", "metabolites"). If the user select "subject", the boxplot will be created for each subject. If the selection is "metabolites", the boxplot will be created for each metabolite.
#'
#' @export
#'
#' @return A ggplot2 object.
#' @author Pol Castellano-Escuder
PomaNormPlot <- function(data, group = c("subjects", "metabolites")){

  if (!(group %in% c("subjects", "metabolites"))) {
    stop(crayon::red(clisymbols::symbol$cross, "Incorrect value for group argument!"))
  }
  if (missing(group)) {
    group <- "subjects"
    warning("group argument is empty! subjects will be used")
  }

  colnames(data)[1:2] <- c("ID", "Group")
  data <- data %>% mutate(ID = paste0(row_number(), "_", ID))

  normtable_metabolites <- dplyr::select(data, -ID) %>%
    reshape2::melt()

  normtable_subjects <- dplyr::select(data, -Group) %>%
    as_tibble() %>%
    gather(var, value, -ID) %>%
    spread(ID, value) %>%
    reshape2::melt()

  normtable_subjects <- normtable_subjects %>%
    mutate(variable = str_replace_all(variable, "._", ""))

  if(group == "subjects"){

    normtable_subjects %>%
      group_by(variable) %>%
      ggplot(aes(variable, value, color = variable)) +
      geom_boxplot() +
      geom_jitter() +
      theme(legend.position = "none") +
      theme_minimal() +
      xlab("") +
      ggtitle("Normalization Plot by Subjects") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "none")

  }
  else{

    normtable_metabolites %>%
      group_by(Group) %>%
      ggplot(aes(variable, value, color = Group)) +
      geom_boxplot() +
      geom_jitter() +
      theme(legend.position = "none") +
      theme_minimal() +
      xlab("") +
      ggtitle("Normalization Plot by Metabolites") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))

  }

}

