
#' Distribution Plot
#'
#' @description PomaDensity() generates a density plot of not normalized and normalized MS data. This plot can help in the comparison between pre and post normalized data and in the "validation" of the normalization process.
#'
#' @param data A MSnSet object. First `pData` column must be the subject group/type.
#' @param group Groupping factor for the plot. Options are "samples" and "features". Option "samples" (default) will create a density plot for each group and option "features" will create a density plot of each variable.
#' @param feature_name A vector with the name/s of feature/s to plot. If it's NULL (default) a density plot of all variables will be created.
#' @param legend_position Character indicating the legend position. Options are "none", "top", "bottom", "left", and "right".
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
#' @importFrom tidyr pivot_longer
#' @importFrom MSnbase pData exprs featureNames
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
  if (missing(group)) {
    message("group argument is empty! samples will be used")
  }
  if (!is.null(feature_name)) {
    if(!isTRUE(all(feature_name %in% MSnbase::featureNames(data)))){
      stop("At least one feature name not found...")
    }
  }
  if(!(legend_position %in% c("none", "top", "bottom", "left", "right"))) {
    stop("Incorrect value for legend_position argument!")
  }

  e <- t(MSnbase::exprs(data))
  target <- MSnbase::pData(data) %>%
    rownames_to_column("ID") %>%
    rename(Group = 2) %>%
    select(ID, Group)
  
  data <- cbind(target, e)
  
  if(group == "samples"){

    if (is.null(feature_name)){
      
      data %>%
        pivot_longer(cols = -c(ID, Group)) %>%
        ggplot(aes(value, fill = Group)) +
        geom_density(alpha = 0.4) +
        xlab("Value") +
        ylab("Density") +
        theme_bw() +
        theme(legend.title = element_blank(),
              legend.position = legend_position) +
        scale_fill_viridis_d()
      
    } else {
      
      data %>%
        pivot_longer(cols = -c(ID, Group)) %>%
        filter(name %in% feature_name) %>%
        ggplot(aes(value, fill = Group)) +
        geom_density(alpha = 0.4) +
        xlab("Value") +
        ylab("Density") +
        theme_bw() +
        theme(legend.title = element_blank(),
              legend.position = legend_position) +
        scale_fill_viridis_d()
      
    }

  } else {

    if (is.null(feature_name)){

      data %>%
        dplyr::select(-ID) %>%
        pivot_longer(cols = -Group) %>%
        ggplot(aes(value, fill = name)) +
        geom_density(alpha = 0.4) +
        theme_bw() +
        xlab("Value") +
        ylab("Density") +
        theme(legend.position = "none") +
        scale_fill_viridis_d()

    } else {
      
      data %>%
        dplyr::select(-ID) %>%
        pivot_longer(cols = -Group) %>%
        filter(name %in% feature_name) %>%
        ggplot(aes(value, fill = name)) + 
        geom_density(alpha = 0.4) +
        theme_bw() +
        xlab("Value") +
        ylab("Density") +
        theme(legend.title = element_blank(),
              legend.position = legend_position) +
        scale_fill_viridis_d()

    }
  }
}

