
#' Density Plots
#'
#' @description `PomaDensity` generates a density plot for samples and features. This function can be used for data exploration (e.g., comparison between pre and post normalized datasets).
#'
#' @param data A `SummarizedExperiment` object.
#' @param x Character. Options are "samples" (to visualize sample density plots) and "features" (to visualize feature density plots). Default is "samples".
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
#' # Sample density plots
#' st000284 %>%
#' PomaNorm() %>% 
#' PomaDensity(theme_params = list(axistext = "y"))
#' 
#' # All feature density plots
#' st000284 %>% 
#' PomaNorm() %>% 
#' PomaDensity(x = "features", theme_params = list(legend_position = "none"))
#'              
#' # Specific feature density plots
#' st000284 %>% 
#' PomaNorm() %>% 
#' PomaDensity(x = "features", 
#'             feature_name = c("ornithine", "orotate"))
PomaDensity <- function(data,
                        x = "samples",
                        feature_name = NULL,
                        theme_params = list(legend_title = FALSE),
                        ...) {

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
    if(length(feature_name) > 10) {
      stop("Maximum 10 features")
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
  
  if (x == "samples"){
    plot_data <- plot_data %>%
      tidyr::pivot_longer(cols = -c(sample_id, group_factor))
    
    if (!is.null(feature_name)){
      plot_data <- plot_data %>% 
        dplyr::filter(name %in% feature_name)
    }
    
    plot_complete <- plot_data %>% 
      ggplot2::ggplot(ggplot2::aes(value)) +
      {if(grouping_factor)ggplot2::geom_density(ggplot2::aes(fill = group_factor), alpha = 0.5)} +
      {if(!grouping_factor)ggplot2::geom_density(alpha = 0.5)} +
      ggplot2::labs(x = "Value", 
                    y = "Density",
                    fill = NULL) +
      do.call(theme_poma, theme_params) +
      scale_fill_poma_d()

  } 
  
  else if (x == "features") {
    plot_data <- plot_data %>%
      dplyr::select(-sample_id) %>%
      tidyr::pivot_longer(cols = -group_factor)
    
    if (!is.null(feature_name)) {
      plot_data <- plot_data %>% 
        dplyr::filter(name %in% feature_name)
    }

    plot_complete <- plot_data %>% 
      ggplot2::ggplot(ggplot2::aes(value, fill = name)) +
      ggplot2::geom_density(alpha = 0.5) +
      ggplot2::labs(x = "Value",
                    y = "Density",
                    fill = NULL) +
      do.call(theme_poma, theme_params) +
      scale_fill_poma_d()
  }
  
  return(plot_complete)
}

