
#' Dimensionality Reduction with UMAP
#'
#' @description Dimension reduction of the data using the Uniform Manifold Approximation and Projection (UMAP) method. See `?uwot::umap()` for more.
#'
#' @param data A SummarizedExperiment object.
#' @param n_neighbors The size of local neighborhood (sample points) used for manifold approximation. 
#' @param n_components The dimension of the space to embed into. Defaults is 2.
#' @param metric Distance measure method to find nearest neighbors. Options are "euclidean", "cosine", "manhattan", "hamming" and "correlation". See `?uwot::umap()`.
#' @param pca If not NULL (default), reduce data to this number of columns using PCA before UMAP. 
#' @param min_dist The effective minimum distance between embedded points.
#' @param spread The effective scale of embedded points.
#' @param hdbscan_minpts Integer; Minimum size of clusters. See `?dbscan::hdbscan()`.
#' @param show_clusters Logical indicating if clusters computed with HDBSCAN method should be plotted or not.
#' @param show_group Logical indicating if the original sample group from target should be plotted or not.
#' @param legend_position Character indicating the legend position. Options are "none", "top", "bottom", "left", and "right".
#' 
#' @export
#'
#' @return A list with results including plots and tables.
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
                     n_neighbors = NULL,
                     n_components = 2,
                     metric = "euclidean",
                     pca = NULL,
                     min_dist = 0.01,
                     spread = 1,
                     hdbscan_minpts = NULL,
                     show_clusters = FALSE,
                     show_group = FALSE,
                     legend_position = "bottom") {
  
  if (missing(data)) {
    stop("data argument is empty!")
  }
  if(!is(data, "SummarizedExperiment")){
    stop("data is not a SummarizedExperiment object. \nSee POMA::PomaSummarizedExperiment or SummarizedExperiment::SummarizedExperiment")
  }
  if(show_group & (!is(SummarizedExperiment::colData(data)[,1], "character") & 
                   !is(SummarizedExperiment::colData(data)[,1], "factor"))){
    stop("show_group = TRUE expects the first column of your target to be a factor or character")
  }
  if(show_clusters & is.null(hdbscan_minpts)){
    hdbscan_minpts <- 5
    warning("HDBSCAN minPts set to 5 samples")
  }
  if(show_clusters & show_group){
    warning("Both show_clusters and show_group are set to TRUE. Choose only one. Showing the clusters by default.")
  }
  
  umap_data <- t(SummarizedExperiment::assay(data))

  if(is.null(n_neighbors)) {
    n_neighbors <- floor(sqrt(nrow(umap_data))) 
  }
  
  umap_res <- uwot::umap(umap_data,
                         n_neighbors = n_neighbors, 
                         n_components = n_components,
                         metric = metric,
                         pca = pca,
                         min_dist = min_dist,
                         spread = spread)
  
  if(!is.null(hdbscan_minpts)) {
    hdbscan_res <- dbscan::hdbscan(umap_res,
                                   minPts = hdbscan_minpts)
    
    umap_result <- dplyr::tibble(sample = rownames(SummarizedExperiment::colData(data)),
                                 as.data.frame(umap_res),
                                 clust = as.factor(hdbscan_res$cluster),
                                 member_prob = hdbscan_res$membership_prob) 
  } else {
    umap_result <- dplyr::tibble(sample = rownames(SummarizedExperiment::colData(data)),
                                 as.data.frame(umap_res))
  }
  
  umap_result %<>%
    dplyr::rename_at(dplyr::vars(dplyr::starts_with("V")), ~ gsub("V", "UMAP", .))
  
  if(show_group) {
    umap_result %<>%
      dplyr::mutate(group = SummarizedExperiment::colData(data)[,1]) %>% 
      dplyr::relocate(group, .after = sample)
  }
    
  umap_plot <- ggplot2::ggplot(umap_result, ggplot2::aes(UMAP1, UMAP2)) +
    {if(!show_group & !show_clusters)ggplot2::geom_point(size = 2)} +
    {if(show_group & !show_clusters)ggplot2::geom_point(ggplot2::aes(color = group), size = 2)} +
    {if(!show_group & show_clusters)ggplot2::geom_point(ggplot2::aes(color = clust), size = 2)} +
    {if(show_group & show_clusters)ggplot2::geom_point(ggplot2::aes(color = clust), size = 2)} +
    ggplot2::labs(x = "UMAP 1",
                  y = "UMAP 2") +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.title = ggplot2::element_blank(),
                   legend.position = legend_position) +
    ggplot2::scale_color_viridis_d(option = "plasma", end = 0.8)
  
  return(list(umap_table = umap_result,
              umap_plot = umap_plot))

}

