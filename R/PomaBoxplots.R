
#' Boxplots and Violin Plots
#'
#' @description `PomaBoxplots` generates boxplots and violin plots for samples and features. This function can be used for data exploration (e.g., comparison between pre and post normalized datasets).
#'
#' @param data A `SummarizedExperiment` object.
#' @param x Character. Options are "samples" (to visualize sample boxplots) and "features" (to visualize feature boxplots). Default is "samples".
#' @param violin Logical. Indicates if violin plots should be displayed instead of boxplots. Default is FALSE.
#' @param feature_name Character vector. Indicates the feature/s to display. Default is NULL (all features will be displayed).
#' @param theme_params List. Indicates `theme_poma` parameters.
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
#' # Sample boxplots
#' st000284 %>%
#' PomaNorm() %>% 
#' PomaBoxplots(theme_params = list(axistext = "y"))
#' 
#' # Sample violin plots
#' st000284 %>%
#' PomaNorm() %>% 
#' PomaBoxplots(violin = TRUE, theme_params = list(axistext = "y"))
#' 
#' # All feature boxplots
#' st000284 %>% 
#' PomaNorm() %>% 
#' PomaBoxplots(x = "features", theme_params = list(axis_x_rotate = TRUE))
#'              
#' # Specific feature boxplots
#' st000284 %>% 
#' PomaNorm() %>% 
#' PomaBoxplots(x = "features", 
#'              feature_name = c("ornithine", "orotate"))
#'              
#' # Specific feature violin plots
#' st000284 %>% 
#' PomaNorm() %>% 
#' PomaBoxplots(x = "features", 
#'              violin = TRUE,
#'              feature_name = c("ornithine", "orotate"))
PomaBoxplots <- function(data,
                         x = "samples",
                         violin = FALSE,
                         feature_name = NULL,
                         theme_params = list(legend_title = FALSE, axis_x_rotate = TRUE)) {

  if(!is(data, "SummarizedExperiment")){
    stop("data is not a SummarizedExperiment object. \nSee POMA::PomaCreateObject or SummarizedExperiment::SummarizedExperiment")
  }
  if (!(x %in% c("samples", "features"))) {
    stop("Incorrect value for x argument")
  }
  if (!is.null(feature_name)) {
    if(!any(feature_name %in% rownames(SummarizedExperiment::assay(data)))) {
      stop("Features not found")
    }
    if(!all(feature_name %in% rownames(SummarizedExperiment::assay(data)))) {
      message(paste0(paste0(feature_name[!feature_name %in% rownames(SummarizedExperiment::assay(data))], collapse = ", "),
                     " not found"))
    }
  }
  
  plot_data <- t(SummarizedExperiment::assay(data))
  
  grouping_factor <- ifelse(ncol(SummarizedExperiment::colData(data)) > 0, 
                            is.factor(SummarizedExperiment::colData(data)[,1]), FALSE)
  
  if (grouping_factor) {
    plot_data <- data.frame(sample_id = rownames(SummarizedExperiment::colData(data)),
                            group_factor = SummarizedExperiment::colData(data)[,1], 
                            plot_data)
  } else {
    if (ncol(SummarizedExperiment::colData(data)) > 0) {
      sample_names <- rownames(SummarizedExperiment::colData(data))
    } else {
      sample_names <- colnames(SummarizedExperiment::assay(data))
    }
      
    if (is.null(sample_names)) {
      sample_names <- paste0("sample_", 1:ncol(SummarizedExperiment::assay(data)))
    }
    
    plot_data <- data.frame(sample_id = sample_names,
                            group_factor = "no_groups",
                            plot_data)
  }

  if (x == "samples") {
    plot_data <- plot_data %>%
      tidyr::pivot_longer(cols = -c(sample_id, group_factor)) %>%
      dplyr::group_by(sample_id) %>% 
      dplyr::mutate(median_rank = median(value, na.rm = TRUE)) %>% 
      dplyr::ungroup() %>% 
      ggplot2::ggplot(ggplot2::aes(reorder(sample_id, median_rank), value))
  }
  
  else if (x == "features") {
    plot_data <- plot_data %>%
      dplyr::select(-sample_id) %>%
      tidyr::pivot_longer(cols = -group_factor)
    
    if (!is.null(feature_name)){
      plot_data <- plot_data %>% 
        dplyr::filter(name %in% feature_name)
    }
    
    plot_data <- plot_data %>% 
      ggplot2::ggplot(ggplot2::aes(name, value))
  }
  
  plot_complete <- plot_data +
    # boxplots
    {if(grouping_factor & !violin)ggplot2::geom_boxplot(ggplot2::aes(color = group_factor, fill = group_factor), alpha = 0.5)} +
    {if(!grouping_factor & !violin)ggplot2::geom_boxplot(alpha = 0.5)} +
    # violin plots
    {if(grouping_factor & violin)ggplot2::geom_violin(ggplot2::aes(color = group_factor, fill = group_factor), alpha = 0.5)} +
    {if(!grouping_factor & violin)ggplot2::geom_violin(alpha = 0.5)} +
    # aesthetics
    ggplot2::labs(x = NULL, 
                  y = "Value",
                  fill = NULL,
                  color = NULL) +
    do.call(theme_poma, theme_params) +
    scale_fill_poma_d() +
    scale_color_poma_d()
  
  return(plot_complete)
}

