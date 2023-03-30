
#' Distribution Plot
#'
#' @description PomaDensity() generates a density plot of not normalized and normalized MS data. This plot can help in the comparison between pre and post normalized data and in the "validation" of the normalization process.
#'
#' @param data A SummarizedExperiment object.
#' @param group Groupping factor for the plot. Options are "samples" and "features". Option "samples" (default) will create a density plot for each group and option "features" will create a density plot of each variable.
#' @param feature_name A vector with the name/s of feature/s to plot. If it's NULL (default) a density plot of all variables will be created.
#' @param theme_params List indicating `theme_poma` parameters.
#' @param palette POMA palette. One of "nature", "simpsons", or "futurama".
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
#' PomaDensity(st000284)
#' 
#' # features
#' PomaDensity(st000284, group = "features")
#' 
#' # concrete features
#' PomaDensity(st000284, group = "features", 
#'             feature_name = c("ornithine", "orotate"))
PomaDensity <- function(data,
                        group = "samples",
                        feature_name = NULL,
                        theme_params = list(),
                        palette = "nature",
                        ...) {

  if (missing(data)) {
    stop("data argument is empty!")
  }
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
    if(length(feature_name) > 10) {
      stop("Maximum 10 features")
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
      tidyr::pivot_longer(cols = -c(ID, Group))
    
    if (!is.null(feature_name)){
      plot_data <- plot_data %>% 
        dplyr::filter(name %in% feature_name)
    }
    
    plot_complete <- plot_data %>% 
      ggplot2::ggplot(ggplot2::aes(value, fill = Group)) +
      ggplot2::geom_density(alpha = 0.5) +
      ggplot2::labs(x = "Value", y = "Density") +
      do.call(theme_poma, theme_params) +
      # scale_fill_poma_d(palette = palette)
      ggplot2::scale_fill_viridis_d(option = "plasma", end = 0.8)

  } else {
    plot_data <- data %>%
      dplyr::select(-ID) %>%
      tidyr::pivot_longer(cols = -Group)
    
    if (!is.null(feature_name)) {
      plot_data <- plot_data %>% 
        dplyr::filter(name %in% feature_name)
      }

    plot_data <- plot_data %>% ggplot2::ggplot(ggplot2::aes(value, fill = name))
    
    plot_complete <- plot_data +
      ggplot2::geom_density(alpha = 0.5) +
      ggplot2::labs(x = "Value",
                    y = "Density") +
      do.call(theme_poma, theme_params)
  }
  
  return(plot_complete)
  
}

