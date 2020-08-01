
#' Classical Boxplots
#'
#' @description PomaBoxplots() generates a boxplot for subjects or features. This boxplot can help in the comparison between pre and post normalized data and in the "validation" of the normalization process.
#'
#' @param data A MSnSet object. First `pData` column must be the subject group/type.
#' @param group Groupping factor for the plot. Options are "samples" and "features". Option "samples" (default) will create a boxplot for each sample and option "features" will create a boxplot of each variable.
#' @param jitter Logical. If it's TRUE (default), the boxplot will show all points.
#' @param feature_name A vector with the name/s of feature/s to plot. If it's NULL (default) a boxplot of all features will be created.
#' @param label_size Numeric indicating the size of x-axis labels.
#'
#' @export
#'
#' @return A ggplot2 object.
#' @author Pol Castellano-Escuder
#'
#' @import ggplot2
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr select group_by filter rename
#' @importFrom magrittr %>%
#' @importFrom reshape2 melt
#' @importFrom crayon red
#' @importFrom clisymbols symbol
#' @importFrom Biobase varLabels pData exprs featureNames
#' 
#' @examples 
#' data("st000284")
#' 
#' # samples
#' PomaBoxplots(st000284)
#' 
#' # features
#' PomaBoxplots(st000284, group = "features")
#'              
#' # concrete features
#' PomaBoxplots(st000284, group = "features", 
#'              feature_name = c("ornithine", "orotate"))
PomaBoxplots <- function(data,
                         group = "samples",
                         jitter = TRUE,
                         feature_name = NULL,
                         label_size = 10){
  
  if(missing(data)) {
    stop(crayon::red(clisymbols::symbol$cross, "data argument is empty!"))
  }
  if(!is(data[1], "MSnSet")){
    stop(paste0(crayon::red(clisymbols::symbol$cross, "data is not a MSnSet object."), 
                " \nSee POMA::PomaMSnSetClass or MSnbase::MSnSet"))
  }
  if (!(group %in% c("samples", "features"))) {
    stop(crayon::red(clisymbols::symbol$cross, "Incorrect value for group argument!"))
  }
  if (!is.null(feature_name)) {
    if(!isTRUE(all(feature_name %in% Biobase::featureNames(data)))){
      stop(crayon::red(clisymbols::symbol$cross, "At least one feature name not found..."))
    }
  }
  
  e <- t(Biobase::exprs(data))
  target <- Biobase::pData(data) %>%
    rownames_to_column("ID") %>%
    rename(Group = 2) %>%
    select(ID, Group)
  
  data <- cbind(target, e)

  if(group == "samples"){
    
    data %>%
      reshape2::melt() %>%
      group_by(ID) %>%
      ggplot(aes(ID, value, color = Group)) +
      geom_boxplot() +
      {if(jitter)geom_jitter(alpha = 0.5, position = position_jitterdodge())} +
      theme_bw() +
      xlab("") +
      ylab("Value") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = label_size))
  }
  
  else {
    
    if(is.null(feature_name)){
      
      data %>%
        dplyr::select(-ID) %>%
        reshape2::melt() %>%
        group_by(Group) %>%
        ggplot(aes(variable, value, color = Group)) +
        geom_boxplot() +
        {if(jitter)geom_jitter(alpha = 0.5, position = position_jitterdodge())} +
        theme_bw() +
        xlab("") +
        ylab("Value") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size = label_size))
      
    } else {
      
      data %>%
        dplyr::select(-ID) %>%
        reshape2::melt() %>%
        group_by(Group) %>%
        filter(variable %in% feature_name) %>%
        ggplot(aes(variable, value, color = Group)) +
        geom_boxplot() +
        {if(jitter)geom_jitter(alpha = 0.5, position = position_jitterdodge())} +
        theme_bw() +
        xlab("") +
        ylab("Value") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size = label_size))

    }
  }
}

