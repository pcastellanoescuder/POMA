
#' Classical Boxplots
#'
#' @description PomaBoxplots() generates a boxplot for subjects or features. This boxplot can help in the comparison between pre and post normalized data and in the "validation" of the normalization process.
#'
#' @param data A SummarizedExperiment object.
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
#' @importFrom magrittr %>%
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
    tibble::rownames_to_column("ID") %>%
    dplyr::rename(Group = 2) %>%
    dplyr::select(ID, Group)
  
  data <- cbind(target, e)

  if(group == "samples"){
    plot_data <- data %>%
      tidyr::pivot_longer(cols = -c(ID, Group)) %>%
      ggplot2::ggplot(ggplot2::aes(ID, value, color = Group))
  }
  
  else {
    if(is.null(feature_name)){
      plot_data <- data %>%
        dplyr::select(-ID) %>%
        tidyr::pivot_longer(cols = -Group) %>%
        ggplot2::ggplot(ggplot2::aes(name, value, color = Group))
      
    } else {
      plot_data <- data %>%
        dplyr::select(-ID) %>%
        tidyr::pivot_longer(cols = -Group) %>%
        dplyr::filter(name %in% feature_name) %>%
        ggplot2::ggplot(ggplot2::aes(name, value, color = Group))
    }
  }
  
  plot_complete <- plot_data +
    ggplot2::geom_boxplot() +
    {if(jitter)ggplot2::geom_jitter(alpha = 0.5, position = ggplot2::position_jitterdodge())} +
    ggplot2::theme_bw() +
    ggplot2::labs(x = "", 
                  y = "Value") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = label_size),
                   legend.title = ggplot2::element_blank(),
                   legend.position = legend_position) +
    ggplot2::scale_colour_viridis_d(begin = 0, end = 0.8)
  
  return(plot_complete)
  
}

