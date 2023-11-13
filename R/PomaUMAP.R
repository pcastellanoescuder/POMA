
#' Dimensionality Reduction with UMAP
#'
#' @description `PomaUMAP` performs a dimension reduction of the data using the Uniform Manifold Approximation and Projection (UMAP) method. See `?uwot::umap()` for more.
#'
#' @param data A `SummarizedExperiment` object.
#' @param n_neighbors Numeric. Indicates the size of local neighborhood (sample points) used for manifold approximation. 
#' @param n_components Numeric. Indicates the dimension of the space to embed into.
#' @param metric Character. Indicates the distance measure method to find nearest neighbors. Options are "euclidean", "cosine", "manhattan", "hamming" and "correlation". See `?uwot::umap()`.
#' @param pca If not NULL (default), reduce data to this number of columns using PCA before UMAP. 
#' @param min_dist Numeric. Indicates the effective minimum distance between embedded points.
#' @param spread Numeric. Indicates the effective scale of embedded points.
#' @param hdbscan_minpts Numeric. Indicates the minimum size of clusters. See `?dbscan::hdbscan()`.
#' @param show_clusters Logical. Indicates if clusters computed with HDBSCAN method should be plotted or not.
#' @param hide_noise Logical. Specifies whether to hide Cluster 0 in the plot. In HDBSCAN, Cluster 0 is typically regarded as "noise."
#' @param theme_params List. Indicates `theme_poma` parameters.
#' 
#' @export
#'
#' @return A `list` with results including plots and tables.
#' @references McInnes, L., Healy, J., & Melville, J. (2018). Umap: Uniform manifold approximation and projection for dimension reduction. arXiv preprint arXiv:1802.03426.
#' @references Campello, R. J., Moulavi, D., & Sander, J. (2013, April). Density-based clustering based on hierarchical density estimates. In Pacific-Asia conference on knowledge discovery and data mining (pp. 160-172). Springer, Berlin, Heidelberg.
#' @author Pol Castellano-Escuder
#'
#' @importFrom magrittr %>% %<>%
#' 
#' @examples 
#' data("st000284")
#' 
#' st000284 %>%
#'   PomaNorm() %>%
#'   PomaUMAP()
PomaUMAP <- function(data,
                     n_neighbors = floor(sqrt(nrow(data))),
                     n_components = 2,
                     metric = "euclidean",
                     pca = NULL,
                     min_dist = 0.01,
                     spread = 1,
                     hdbscan_minpts = floor(nrow(data) * 0.05),
                     show_clusters = TRUE,
                     hide_noise = TRUE,
                     theme_params = list(legend_title = TRUE, legend_position = "bottom")) {
  
  if(!is(data, "SummarizedExperiment")){
    stop("data is not a SummarizedExperiment object. \nSee POMA::PomaCreateObject or SummarizedExperiment::SummarizedExperiment")
  }
  
  to_umap <- t(SummarizedExperiment::assay(data))
  
  umap_res <- uwot::umap(to_umap,
                         n_neighbors = n_neighbors, 
                         n_components = n_components,
                         metric = metric,
                         pca = pca,
                         min_dist = min_dist,
                         spread = spread)
  
  hdbscan_res <- dbscan::hdbscan(umap_res, minPts = hdbscan_minpts)
  
  umap_clusters <- dplyr::tibble(sample = rownames(SummarizedExperiment::colData(data)),
                                 as.data.frame(umap_res),
                                 clust = as.factor(hdbscan_res$cluster),
                                 member_prob = hdbscan_res$membership_prob) 

  umap_clusters %<>%
    dplyr::rename_at(dplyr::vars(dplyr::starts_with("V")), ~ gsub("V", "UMAP", .))
  
  if (hide_noise) {
    plot_data <- umap_clusters %>% 
      dplyr::filter(clust != 0)
  } else {
    plot_data <- umap_clusters
  }
  
  umap_plot <- ggplot2::ggplot(plot_data, ggplot2::aes(UMAP1, UMAP2)) +
    {if(!show_clusters)ggplot2::geom_point(pch = 21, size = 3, alpha = 0.8)} +
    {if(show_clusters)ggplot2::geom_point(ggplot2::aes(fill = clust), pch = 21, size = 3, alpha = 0.8)} +
    ggplot2::labs(x = "UMAP 1",
                  y = "UMAP 2",
                  fill = "Cluster") +
    do.call(theme_poma, theme_params) +
    POMA::scale_fill_poma_d()

  return(list(umap_embeddings = umap_clusters,
              umap_plot = umap_plot))
}

