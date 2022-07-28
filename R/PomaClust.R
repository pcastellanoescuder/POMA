
#' Cluster Analysis
#'
#' @description This function performs a classical multidimensional scaling (MDS) using all features in the data and computes a cluster analysis for `k` clusters. Then, the calculated clusters will be represented on a MDS plot.
#' 
#' @param data A SummarizedExperiment object.
#' @param method Distance measure method to perform MDS. Options are "euclidean", "maximum", "manhattan", "canberra" and "minkowski". See `?dist()`.
#' @param k Number of clusters (default is `NA`). The optimum number of clusters will be used by default.
#' @param k_max Number of clusters among which the optimal one will be selected.
#' @param show_clusters Logical indicating if clusters should be plotted or not. If this parameter is set to `FALSE` the resultant plot will be a classical 2-dimension MDS plot.
#' @param labels Logical indicating if sample names should be plotted or not.
#' @param show_group Logical indicating if the original sample group from pData should be plotted instead of sample ID or not. Only works if labels is set to `TRUE`.
#' 
#' @export
#'
#' @return A list with the results.
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
                      labels = FALSE,
                      show_group = FALSE){
  
  if (missing(data)) {
    stop("data argument is empty!")
  }
  if(!is(data[1], "SummarizedExperiment")){
    stop("data is not a SummarizedExperiment object. \nSee POMA::PomaSummarizedExperiment or SummarizedExperiment::SummarizedExperiment")
  }
  if (!(method %in% c("euclidean", "maximum", "manhattan", "canberra", "minkowski"))) {
    stop("Incorrect value for method argument!")
  }
  
  e <- t(SummarizedExperiment::assay(data))
  target <- SummarizedExperiment::colData(data)
  
  ## Optimal number of clusters
  wss <- data.frame(wss = sapply(1:k_max, function(k){stats::kmeans(e, k)$tot.withinss})) %>%
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
    ggplot2::geom_point(ggplot2::aes(x = k[elbowi], y = wss[elbowi]), color = "red", size = 2) +
    ggplot2::theme_bw()
  
  ## MDS and kmeans
  if(is.na(k)){
    k <- elbowi
  }
  
  clusters <- stats::kmeans(e, centers = k)
  
  mds <- e %>% 
    stats::dist(method = method) %>% 
    stats::cmdscale() %>% 
    as.data.frame() %>% 
    dplyr::mutate(sample = rownames(target),
                  group = target[,1],
                  clust = as.factor(clusters$cluster)) %>%
    dplyr::select(sample, group, clust, V1, V2) %>% 
    dplyr::rename(Dim1 = V1, Dim2 = V2) %>% 
    dplyr::as_tibble()
  
  if(!show_clusters){
    
    mds_plot <- ggplot2::ggplot(mds, ggplot2::aes(x = Dim1, y = Dim2)) +
      {if(!labels)ggplot2::geom_point(size = 3, alpha = 0.8)} +
      ggplot2::labs(x = "Dimension 1",
                    y = "Dimension 2") +
      ggplot2::theme_bw() +
      {if(labels & show_group)ggplot2::geom_text(ggplot2::aes(label = group))} +
      {if(labels & !show_group)ggplot2::geom_text(ggplot2::aes(label = sample))}
    
  } else {
    
    mds_plot <- ggplot2::ggplot(mds, ggplot2::aes(x = Dim1, y = Dim2, color = clust)) +
      {if(!labels)ggplot2::geom_point(size = 3, alpha = 0.8)} +
      ggplot2::labs(x = "Dimension 1",
                    y = "Dimension 2") +
      ggplot2::theme_bw() +
      {if(labels & show_group)ggplot2::geom_text(ggplot2::aes(label = group), show.legend = FALSE)} +
      {if(labels & !show_group)ggplot2::geom_text(ggplot2::aes(label = sample), show.legend = FALSE)} +
      ggplot2::scale_color_viridis_d(begin = 0, end = 0.8)
  }
  
  return(list(mds_values = mds, 
              optimum_cluster_num = elbowi, 
              optimum_cluster_plot = optimum_clusters, 
              mds_plot = mds_plot))
  
}

