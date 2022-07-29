
#' Dimensionality Reduction with UMAP
#'
#' @description Carry out dimensionality reduction of a dataset using the Uniform Manifold Approximation and Projection (UMAP) method. 
#'
#' @param data A SummarizedExperiment object.
#'
#' @export
#'
#' @return A list with results.
#' @references McInnes, L., Healy, J., & Melville, J. (2018). Umap: Uniform manifold approximation and projection for dimension reduction. arXiv preprint arXiv:1802.03426.
#' @author Pol Castellano-Escuder
#'
#' @importFrom magrittr %>%
#' 
#' @examples 
#' data("st000284")
#' 
#' st000284 %>%
#'   PomaNorm() %>%
#'   PomaUMAP()
PomaUMAP <- function(data) {
  
  if (missing(data)) {
    stop("data argument is empty!")
  }
  if(!is(data, "SummarizedExperiment")){
    stop("data is not a SummarizedExperiment object. \nSee POMA::PomaSummarizedExperiment or SummarizedExperiment::SummarizedExperiment")
  }

  e <- SummarizedExperiment::assay(data)
  
}

