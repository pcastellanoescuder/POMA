
#' Analyse and Remove Statistical Outliers
#'
#' @description `PomaOutliers` analyses and removes statistical outliers from the data.
#'
#' @param data A `SummarizedExperiment` object.
#' @param method Character. Indicates the distance measure method to perform MDS.
#' @param type Character. Indicates the type of outlier analysis to perform. Options are "median" (default) and "centroid". See `vegan::betadisper`.
#' @param coef Numeric. Indicates the outlier coefficient. Lower values are more sensitive to outliers while higher values are less restrictive about outliers. 
#' @param labels Logical. Indicates if sample names should to be plotted.
#'
#' @export
#'
#' @return A `list` with the results.
#' @author Pol Castellano-Escuder
#'
#' @importFrom magrittr %>%
#' 
#' @examples 
#' data("st000336")
#' 
#' # clean outliers
#' st000336 %>% 
#'   PomaImpute() %>%
#'   PomaNorm() %>%
#'   PomaOutliers()
PomaOutliers <- function(data,
                         method = "euclidean",
                         type = "median",
                         coef = 2,
                         labels = FALSE) {
  
  if(!is(data, "SummarizedExperiment")) {
    stop("data is not a SummarizedExperiment object. \nSee POMA::PomaCreateObject or SummarizedExperiment::SummarizedExperiment")
  }
  if (!(method %in% c("euclidean", "maximum", "manhattan", "canberra", "minkowski"))) {
    stop("Incorrect value for method argument")
  }
  if (!(type %in% c("median", "centroid"))) {
    stop("Incorrect value for type argument")
  }
  
  if (ncol(SummarizedExperiment::colData(data)) > 0) {
    group_factor <- SummarizedExperiment::colData(data)[,1]
  } else {
    group_factor <- rep("All Samples", ncol(SummarizedExperiment::assay(data)))
  }
  
  to_outliers <- t(SummarizedExperiment::assay(data))
  
  dd <- stats::dist(to_outliers, method = method)
  distances <- vegan::betadisper(dd, group_factor, type = type, bias.adjust = FALSE, sqrt.dist = FALSE, add = FALSE)
  detect_outliers <- data.frame(distances = distances$distances, 
                                groups = distances$group, 
                                sample = rownames(SummarizedExperiment::colData(data))) %>% 
    dplyr::group_by(groups) %>% 
    dplyr::mutate(limit = quantile(distances, 0.75) + coef * IQR(distances),
                  out = ifelse(distances > limit, 1, 0))
  
  final_outliers <- detect_outliers %>%
    dplyr::filter(out == 1) %>%
    dplyr::select(sample, groups, distance_to_centroid = distances, limit_distance = limit) %>%
    dplyr::as_tibble()
  
  # Analyse
  vectors <- data.frame(distances$vectors)
  centroids <- data.frame(distances$centroids)
  
  total_outliers <- vectors %>%
    tibble::rownames_to_column("sample") %>%
    dplyr::mutate(group = group_factor)
  
  find_hull <- function(x){x[chull(x$PCoA1, x$PCoA2) ,]}
  
  hulls <- total_outliers %>% 
    dplyr::select(-sample) %>%
    dplyr::group_by(group) %>% 
    dplyr::do(find_hull(.)) %>% 
    dplyr::ungroup()
  
  polygon_plot <- ggplot2::ggplot(total_outliers, ggplot2::aes(x = PCoA1, y = PCoA2)) +
    ggplot2::geom_polygon(data = hulls, alpha = 0.4, ggplot2::aes(fill = group)) +
    {if(!labels)ggplot2::geom_point(ggplot2::aes(fill = group), size = 3, alpha = 0.8, pch = 21)} +
    ggplot2::geom_label(data = centroids, ggplot2::aes(x = PCoA1, y = PCoA2, color = rownames(centroids), 
                                                       label = rownames(centroids)), show.legend = FALSE) +
    {if(labels)ggplot2::geom_text(ggplot2::aes(label = sample))} +
    ggplot2::labs(x = "Coordinate 1",
                  y = "Coordinate 2") +
    theme_poma(legend_title = FALSE) +
    scale_fill_poma_d() +
    scale_color_poma_d()
  
  distance_boxplot <- ggplot2::ggplot(detect_outliers, ggplot2::aes(groups, distances, fill = groups)) +
    ggplot2::geom_boxplot(coef = coef, alpha = 0.8) +
    ggplot2::labs(x = NULL,
                  y = "Distance to group centroid") + 
    {if(labels)ggrepel::geom_label_repel(data = detect_outliers[detect_outliers$out == 1,], 
                                         ggplot2::aes(label = sample), alpha = 0.7, na.rm = TRUE, size = 4, show.legend = FALSE)} +
    theme_poma(legend_title = FALSE) +
    scale_fill_poma_d()
  
  # Remove outliers
  metadata <- SummarizedExperiment::colData(data) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("sample") %>%
    dplyr::filter(!(sample %in% final_outliers$sample))
  
  to_outliers_clean <- to_outliers %>%
    as.data.frame() %>%
    dplyr::filter(!(rownames(.) %in% final_outliers$sample))
  
  data_clean <- PomaCreateObject(features = to_outliers_clean, metadata = metadata)
  
  if (validObject(data_clean))
    return(list(polygon_plot = polygon_plot,
                distance_boxplot = distance_boxplot,
                outliers = final_outliers,
                data = data_clean)
           )
}

