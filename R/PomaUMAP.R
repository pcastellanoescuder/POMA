
#' Dimensionality Reduction with UMAP
#'
#' @description Dimension reduction of the data using the Uniform Manifold Approximation and Projection (UMAP) method. 
#'
#' @param data A SummarizedExperiment object.
#' @param n_neighbors xxx
#' @param n_components xxx
#' @param metric xxx
#' @param pca xxx
#' @param min_dist xxx
#' @param spread xxx
#' @param hdbscan_minpts xxx
#' 
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
                     hdbscan_minpts = NULL) {
  
  if (missing(data)) {
    stop("data argument is empty!")
  }
  if(!is(data, "SummarizedExperiment")){
    stop("data is not a SummarizedExperiment object. \nSee POMA::PomaSummarizedExperiment or SummarizedExperiment::SummarizedExperiment")
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
    
  umap_plot <- ggplot2::ggplot(umap_result, ggplot2::aes(UMAP1, UMAP2)) +
    ggplot2::geom_point(size = 2) +
    ggplot2::labs(x = "UMAP 1",
                  y = "UMAP 2") +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "none")
  
  return(list(umap_table = umap_result,
              umap_plot = umap_plot))

}

