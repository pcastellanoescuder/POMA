
#' Cluster Analysis
#'
#' @description This function performs a classical multidimensional scaling (MDS) using all features in the data and computes a cluster analysis for `k` clusters. Then, the calculated clusters will be represented on a MDS plot.
#' 
#' @param data A MSnSet object. First `pData` column must be the subject group/type.
#' @param method Distance measure method to perform MDS. Options are "euclidean", "maximum", "manhattan", "canberra", "binary" and "minkowski". See `?dist()`.
#' @param k Number of clusters.
#' @param show_clusters Logical indicating if clusters should be plotted or not. If this parameter is set to `FALSE` the resultant plot will be a classical 2-dimension MDS plot.
#' @param show_labels Logical indicating if sample names should be plotted or not.
#' 
#' @export
#'
#' @return A list with the results.
#' @author Pol Castellano-Escuder
#'
#' @import ggplot2
#' @importFrom crayon red
#' @importFrom clisymbols symbol
#' @importFrom Biobase exprs pData
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate
PomaClust <- function(data,
                      method = "euclidean",
                      k = 3,
                      show_clusters = TRUE,
                      show_labels = FALSE){
  
  if (missing(data)) {
    stop(crayon::red(clisymbols::symbol$cross, "data argument is empty!"))
  }
  if(!(class(data) == "MSnSet")){
    stop(paste0(crayon::red(clisymbols::symbol$cross, "data is not a MSnSet object."), 
                " \nSee POMA::PomaMSnSetClass or MSnbase::MSnSet"))
  }
  if (!(method %in% c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski"))) {
    stop(crayon::red(clisymbols::symbol$cross, "Incorrect value for method argument!"))
  }
  
  e <- t(Biobase::exprs(data))
  target <- Biobase::pData(data)
  
  ## MDS and kmeans
  
  clusters <- kmeans(e, centers = k)
  
  mds <- e %>% 
    dist(method = method) %>% 
    cmdscale() %>% 
    as.data.frame() %>% 
    mutate(sample = rownames(target),
           cluster = as.factor(clusters$cluster))
  colnames(mds) <- c("Dim1", "Dim2", "sample", "cluster")
    
  ## plot
  
  if(!isTRUE(show_clusters)){
    
    mds_plot <- ggplot(mds, aes(x = Dim1, y = Dim2, label = sample)) +
      geom_point(size = 3, alpha = 0.5) +
      xlab("Dimension 1") +
      ylab("Dimension 2") +
      theme_bw() +
      {if(show_labels)geom_text()}
    
  } else{
    
    mds_plot <- ggplot(mds, aes(x = Dim1, y = Dim2, label = sample, color = cluster, shape = cluster)) +
      geom_point(size = 3, alpha = 0.5) +
      xlab("Dimension 1") +
      ylab("Dimension 2") +
      theme_bw() +
      {if(show_labels)geom_text(show.legend = FALSE)}
  }
  
  return(list(mds_values = mds, mds_plot = mds_plot))
  
}

