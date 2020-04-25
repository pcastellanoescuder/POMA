
#' Classical Boxplots
#'
#' @description PomaBoxplots() generates a boxplot for subjects or features. This boxplot can help in the comparison between pre and post normalized data and in the "validation" of the normalization process.
#'
#' @param data A MSnSet object. First `pData` column must be the subject group/type.
#' @param group Groupping factor for the plot. Options are "samples" and "features". Option "samples" (default) will create a boxplot for each sample and option "features" will create a boxplot of each variable.
#' @param jitter Logical. If it's TRUE (default), the boxplot will show all points.
#' @param feature_name A vector with the name/s of feature/s to plot. If it's NULL (default) a boxplot of all features will be created.
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
#' @importFrom Biobase varLabels pData exprs
PomaBoxplots <- function(data,
                         group = "samples",
                         jitter = TRUE,
                         feature_name = NULL){
  
  if (missing(data)) {
    stop(crayon::red(clisymbols::symbol$cross, "data argument is empty!"))
  }
  if(!(class(data) == "MSnSet")){
    stop(paste0(crayon::red(clisymbols::symbol$cross, "data is not a MSnSet object."), 
                " \nSee POMA::PomaMSnSetClass or MSnbase::MSnSet"))
  }
  if (!(group %in% c("samples", "features"))) {
    stop(crayon::red(clisymbols::symbol$cross, "Incorrect value for group argument!"))
  }
  if (missing(group)) {
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
      xlab("") +
      ylab("Value") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  }
  
  else {
    
    if (is.null(feature_name)){
      
      data %>%
        dplyr::select(-ID) %>%
        reshape2::melt() %>%
        group_by(Group) %>%
        ggplot(aes(variable, value, color = Group)) +
        geom_boxplot() +
        {if(jitter)geom_jitter()} +
        theme_bw() +
        xlab("") +
        ylab("Value") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
      
    } else {
      
      data %>%
        dplyr::select(-ID) %>%
        reshape2::melt() %>%
        group_by(Group) %>%
        filter(variable %in% feature_name) %>%
        ggplot(aes(variable, value, color = Group)) +
        geom_boxplot() +
        {if(jitter)geom_jitter()} +
        theme_bw() +
        xlab("") +
        ylab("Value") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))

    }
  }
}

