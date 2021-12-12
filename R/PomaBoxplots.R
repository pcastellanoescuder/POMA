
#' Classical Boxplots
#'
#' @description PomaBoxplots() generates a boxplot for subjects or features. This boxplot can help in the comparison between pre and post normalized data and in the "validation" of the normalization process.
#'
#' @param data A SummarizedExperiment object. First `colData` column must be the subject group/type.
#' @param group Groupping factor for the plot. Options are "samples" and "features". Option "samples" (default) will create a boxplot for each sample and option "features" will create a boxplot of each variable.
#' @param jitter Logical. If it's TRUE (default), the boxplot will show all points.
#' @param feature_name A vector with the name/s of feature/s to plot. If it's NULL (default) a boxplot of all features will be created.
#' @param label_size Numeric indicating the size of x-axis labels.
#' @param legend_position Character indicating the legend position. Options are "none", "top", "bottom", "left", and "right".
#'
#' @export
#'
#' @return A ggplot2 object.
#' @author Pol Castellano-Escuder
#'
#' @import ggplot2
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr select filter rename
#' @importFrom tidyr pivot_longer
#' @importFrom magrittr %>%
#' @importFrom SummarizedExperiment assay colData
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
                         jitter = FALSE,
                         feature_name = NULL,
                         label_size = 10,
                         legend_position = "bottom"){
  
  if (missing(data)) {
    stop("data argument is empty!")
  }
  if(!is(data[1], "SummarizedExperiment")){
    stop("data is not a SummarizedExperiment object. \nSee POMA::PomaSummarizedExperiment or SummarizedExperiment::SummarizedExperiment")
  }
  if (!(group %in% c("samples", "features"))) {
    stop("Incorrect value for group argument!")
  }
  if (!is.null(feature_name)) {
    if(!isTRUE(all(feature_name %in% rownames(SummarizedExperiment::assay(data))))){
      stop("At least one feature name not found...")
    }
  }
  if(!(legend_position %in% c("none", "top", "bottom", "left", "right"))) {
    stop("Incorrect value for legend_position argument!")
  }
  
  e <- t(SummarizedExperiment::assay(data))
  target <- SummarizedExperiment::colData(data) %>%
    as.data.frame() %>% 
    rownames_to_column("ID") %>%
    rename(Group = 2) %>%
    select(ID, Group)
  
  data <- cbind(target, e)

  if(group == "samples"){
    
    data %>%
      pivot_longer(cols = -c(ID, Group)) %>%
      ggplot(aes(ID, value, color = Group)) +
      geom_boxplot() +
      {if(jitter)geom_jitter(alpha = 0.5, position = position_jitterdodge())} +
      theme_bw() +
      xlab("") +
      ylab("Value") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = label_size),
            legend.title = element_blank(),
            legend.position = legend_position) +
      scale_colour_viridis_d()
  }
  
  else {
    
    if(is.null(feature_name)){
      
      data %>%
        dplyr::select(-ID) %>%
        pivot_longer(cols = -Group) %>%
        ggplot(aes(name, value, color = Group)) +
        geom_boxplot() +
        {if(jitter)geom_jitter(alpha = 0.5, position = position_jitterdodge())} +
        theme_bw() +
        xlab("") +
        ylab("Value") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size = label_size),
              legend.title = element_blank(),
              legend.position = legend_position) +
        scale_colour_viridis_d()
      
    } else {
      
      data %>%
        dplyr::select(-ID) %>%
        pivot_longer(cols = -Group) %>%
        filter(name %in% feature_name) %>%
        ggplot(aes(name, value, color = Group)) +
        geom_boxplot() +
        {if(jitter)geom_jitter(alpha = 0.5, position = position_jitterdodge())} +
        theme_bw() +
        xlab("") +
        ylab("Value") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size = label_size),
              legend.title = element_blank(),
              legend.position = legend_position) +
        scale_colour_viridis_d()

    }
  }
}

