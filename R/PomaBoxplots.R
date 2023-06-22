
#' Classical Boxplots
#'
#' @description PomaBoxplots() generates a boxplot for subjects or features. This boxplot can help in the comparison between pre and post normalized data and in the "validation" of the normalization process.
#'
#' @param data A SummarizedExperiment object.
#' @param group Grouping factor for the plot. Options are "samples" and "features". Option "samples" (default) will create a boxplot for each sample and option "features" will create a boxplot of each variable.
#' @param jitter Logical. If it's TRUE (default), the boxplot will show all points.
#' @param feature_name A vector with the name/s of feature/s to plot. If it's NULL (default) a boxplot of all features will be created.
#' @param theme_params List indicating `theme_poma` parameters.
#' @param palette POMA palette. One of "nature", "simpsons", or "futurama".
#'
#' @export
#'
#' @return A `ggplot` object.
#' @author Pol Castellano-Escuder
#'
#' @importFrom magrittr %>%
#' 
#' @examples 
#' data("st000284")
#' 
#' # samples
#' st000284 %>%
#' PomaNorm() %>% 
#' PomaBoxplots(theme_params = list(axistext = "y"))
#' 
#' # features
#' st000284 %>% 
#' PomaNorm() %>% 
#' PomaBoxplots(group = "features", theme_params = list(axis_x_rotate = TRUE))
#'              
#' # concrete features
#' PomaBoxplots(st000284, group = "features", 
#'              feature_name = c("ornithine", "orotate"))
PomaBoxplots <- function(data,
                         group = "samples",
                         jitter = FALSE,
                         feature_name = NULL,
                         theme_params = list(legend_title = FALSE, axis_x_rotate = TRUE),
                         palette = "nature",
                         ...) {

  if(!is(data, "SummarizedExperiment")){
    stop("data is not a SummarizedExperiment object. \nSee POMA::PomaSummarizedExperiment or SummarizedExperiment::SummarizedExperiment")
  }
  if (!(group %in% c("samples", "features"))) {
    stop("Incorrect value for group argument!")
  }
  if (!is.null(feature_name)) {
    if(!any(feature_name %in% rownames(SummarizedExperiment::assay(data)))) {
      stop("None of the specified features found")
    }
    if(!all(feature_name %in% rownames(SummarizedExperiment::assay(data)))){
      warning(paste0("Feature/s ",
                     paste0(feature_name[!feature_name %in% rownames(SummarizedExperiment::assay(data))], collapse = ", "),
                     " not found"))
    }
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
      ggplot2::ggplot(ggplot2::aes(ID, value, color = Group, fill = Group))
  }
  
  else {
    plot_data <- data %>%
      dplyr::select(-ID) %>%
      tidyr::pivot_longer(cols = -Group)
    
    if(!is.null(feature_name)){
      plot_data <- plot_data %>% 
        dplyr::filter(name %in% feature_name)
    }
    
    plot_data <- plot_data %>% 
      ggplot2::ggplot(ggplot2::aes(name, value, color = Group, fill = Group))
    
  }
  
  plot_complete <- plot_data +
    ggplot2::geom_boxplot(alpha = 0.5) +
    {if(jitter)ggplot2::geom_jitter(alpha = 0.5, position = ggplot2::position_jitterdodge())} +
    ggplot2::labs(x = NULL, y = "Value") +
    do.call(theme_poma, theme_params) +
    # scale_color_poma_d(palette = palette) +
    # scale_fill_poma_d(palette = palette)
    ggplot2::scale_fill_viridis_d(option = "plasma", end = 0.8) +
    ggplot2::scale_color_viridis_d(option = "plasma", end = 0.8)
  
  return(plot_complete)
  
}

