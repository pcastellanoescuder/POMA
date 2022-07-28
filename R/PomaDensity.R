
#' Distribution Plot
#'
#' @description PomaDensity() generates a density plot of not normalized and normalized MS data. This plot can help in the comparison between pre and post normalized data and in the "validation" of the normalization process.
#'
#' @param data A SummarizedExperiment object.
#' @param group Groupping factor for the plot. Options are "samples" and "features". Option "samples" (default) will create a density plot for each group and option "features" will create a density plot of each variable.
#' @param feature_name A vector with the name/s of feature/s to plot. If it's NULL (default) a density plot of all variables will be created.
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
    if (is.null(feature_name)){
      plot_data <- data %>%
        tidyr::pivot_longer(cols = -c(ID, Group)) %>%
        ggplot2::ggplot(ggplot2::aes(value, fill = Group))
      
    } else {
      plot_data <- data %>%
        tidyr::pivot_longer(cols = -c(ID, Group)) %>%
        dplyr::filter(name %in% feature_name) %>%
        ggplot2::ggplot(ggplot2::aes(value, fill = Group))
      
    }

  } else {
    if (is.null(feature_name)){
      plot_data <- data %>%
        dplyr::select(-ID) %>%
        tidyr::pivot_longer(cols = -Group) %>%
        ggplot2::ggplot(ggplot2::aes(value, fill = name))

    } else {
      plot_data <- data %>%
        dplyr::select(-ID) %>%
        tidyr::pivot_longer(cols = -Group) %>%
        dplyr::filter(name %in% feature_name) %>%
        ggplot2::ggplot(ggplot2::aes(value, fill = name))

    }
  }
  
  plot_complete <- plot_data +
    ggplot2::geom_density(alpha = 0.4) +
    ggplot2::theme_bw() +
    ggplot2::labs(x = "Value",
                  y = "Density") +
    ggplot2::theme(legend.title = ggplot2::element_blank(),
                   legend.position = legend_position) +
    ggplot2::scale_fill_viridis_d(begin = 0, end = 0.8)
  
  return(plot_complete)
  
}

