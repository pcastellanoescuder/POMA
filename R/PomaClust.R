
#' Cluster Analysis
#'
#' @description `PomaClust` performs a k-means clustering and plots the results in a classical multidimensional scaling (MDS) plot.
#' 
#' @param data A `SummarizedExperiment` object.
#' @param method Character. Indicates the distance method to perform MDS. Options are "euclidean", "maximum", "manhattan", "canberra" and "minkowski". See `?dist()`.
#' @param k Numeric. Indicates the number of clusters (default is `NA`). The optimal number of clusters will be used by default.
#' @param k_max Numeric. Indicates the number of clusters among which the optimal one will be selected.
#' @param show_clusters Logical. Indicates if clusters should be plotted or not. 
#' @param labels Logical. Indicates if sample names should be plotted or not.
#' 
#' @export
#'
#' @return A `list` with results including plots and tables.
#' @author Pol Castellano-Escuder
#'
#' @importFrom magrittr %>%
#' 
#' @examples 
#' data("st000284")
#' 
#' PomaClust(st000284)
PomaClust <- function(data,
                      method = "euclidean",
                      k = NA,
                      k_max = 15,
                      show_clusters = TRUE,
                      labels = FALSE) {
  
  if (!is(data, "SummarizedExperiment")){
    stop("data is not a SummarizedExperiment object. \nSee POMA::PomaCreateObject or SummarizedExperiment::SummarizedExperiment")
  }
  if (!(method %in% c("euclidean", "maximum", "manhattan", "canberra", "minkowski"))) {
    stop("Incorrect value for method argument")
  }
  
  to_clust <- scale(t(SummarizedExperiment::assay(data)))
  target <- SummarizedExperiment::colData(data)
  
  ## Optimal number of clusters
  wss <- data.frame(wss = sapply(1:k_max, function(k){stats::kmeans(to_clust, k)$tot.withinss})) %>%
    dplyr::mutate(k = 1:k_max)
  
  i1 <- which.min(wss$k)
  i2 <- which.max(wss$k)
  slope <- (wss$wss[i2] - wss$wss[i1]) / (wss$k[i2] - wss$k[i1])
  int <- wss$wss[i1] - slope*wss$k[i1]
  
  perpslope <- -1/slope
  perpint <- wss$wss - perpslope*wss$k
  
  xcross <- (int - perpint) / (perpslope - slope)
  ycross <- slope*xcross + int
  
  dists <- sqrt((wss$k - xcross)^2 + (wss$wss - ycross)^2)
  elbowi <- which.max(dists)
  
  optimum_clusters <- ggplot2::ggplot(wss, ggplot2::aes(as.factor(k), wss, group = 1)) +
    ggplot2::geom_point() +
    ggplot2::geom_line() +
    ggplot2::labs(x = "Number of clusters",
                  y = "Within-cluster sum of squares") +
    ggplot2::geom_point(ggplot2::aes(x = k[elbowi], y = wss[elbowi]), color = "red", size = 3) +
    POMA::theme_poma()
  
  ## MDS and kmeans
  if (is.na(k)) {k <- elbowi}
  
  clusters <- stats::kmeans(to_clust, centers = k)
  
  mds_data <- to_clust %>% 
    stats::dist(method = method) %>% 
    stats::cmdscale() %>% 
    as.data.frame() %>% 
    dplyr::mutate(sample = rownames(target),
                  clust = as.factor(clusters$cluster)) %>%
    dplyr::select(sample, clust, Dim1 = V1, Dim2 = V2) %>% 
    dplyr::as_tibble()
  
  if (!show_clusters) {
    
    mds_plot <- ggplot2::ggplot(mds_data, ggplot2::aes(x = Dim1, y = Dim2)) +
      {if(!labels)ggplot2::geom_point(pch = 21, size = 3, alpha = 0.6)} +
      ggplot2::labs(x = "Dimension 1",
                    y = "Dimension 2") +
      {if(labels)ggplot2::geom_text(ggplot2::aes(label = sample))} +
      POMA::theme_poma()
    
  } else {
    
    small <- nrow(to_clust) < 500
      
    mds_plot <- ggplot2::ggplot(mds_data, ggplot2::aes(x = Dim1, y = Dim2, fill = clust)) +
      ggplot2::geom_point(pch = 21, size = 3, alpha = 0.6) +
      ggplot2::labs(x = "Dimension 1",
                    y = "Dimension 2",
                    fill = "Cluster") +
      {if(labels & !small)ggplot2::geom_text(ggplot2::aes(label = sample), show.legend = FALSE)} +
      {if(labels & small)ggrepel::geom_text_repel(ggplot2::aes(label = sample), show.legend = FALSE)} +
      POMA::theme_poma() +
      POMA::scale_fill_poma_d()
  }
  
  return(list(mds_coordinates = mds_data, 
              mds_plot = mds_plot,
              optimal_clusters_number = elbowi, 
              optimal_clusters_plot = optimum_clusters))
}

