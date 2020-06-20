
#' Cluster Analysis
#'
#' @description This function performs a classical multidimensional scaling (MDS) using all features in the data and computes a cluster analysis for `k` clusters. Then, the calculated clusters will be represented on a MDS plot.
#' 
#' @param data A MSnSet object. First `pData` column must be the subject group/type.
#' @param method Distance measure method to perform MDS. Options are "euclidean", "maximum", "manhattan", "canberra" and "minkowski". See `?dist()`.
#' @param k Number of clusters.
#' @param show_clusters Logical indicating if clusters should be plotted or not. If this parameter is set to `FALSE` the resultant plot will be a classical 2-dimension MDS plot.
#' @param labels Logical indicating if sample names should be plotted or not.
#' @param show_group Logical indicating if the original sample group from pData should be plotted instead of sample ID or not. Only works if labels is set to `TRUE`.
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
#' @importFrom dplyr mutate select
PomaClust <- function(data,
                      method = "euclidean",
                      k = 3,
                      show_clusters = TRUE,
                      labels = FALSE,
                      show_group = FALSE){
  
  if (missing(data)) {
    stop(crayon::red(clisymbols::symbol$cross, "data argument is empty!"))
  }
  if(!(class(data) == "MSnSet")){
    stop(paste0(crayon::red(clisymbols::symbol$cross, "data is not a MSnSet object."), 
                " \nSee POMA::PomaMSnSetClass or MSnbase::MSnSet"))
  }
  if (!(method %in% c("euclidean", "maximum", "manhattan", "canberra", "minkowski"))) {
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
           group = target[,1],
           clust = as.factor(clusters$cluster)) %>%
    dplyr::select(sample, group, clust, V1, V2)
  colnames(mds) <- c("sample", "group", "clust", "Dim1", "Dim2")
    
  ## plot
  
  if(!isTRUE(show_clusters)){
    
    mds_plot <- ggplot(mds, aes(x = Dim1, y = Dim2)) +
      {if(!labels)geom_point(size = 3, alpha = 0.5)} +
      xlab("Dimension 1") +
      ylab("Dimension 2") +
      theme_bw() +
      {if(labels & show_group)geom_text(aes(label = group))} +
      {if(labels & !show_group)geom_text(aes(label = sample))}
    
  } else{
    
    mds_plot <- ggplot(mds, aes(x = Dim1, y = Dim2, color = clust, shape = clust)) +
      {if(!labels)geom_point(size = 3, alpha = 0.5)} +
      xlab("Dimension 1") +
      ylab("Dimension 2") +
      theme_bw() +
      {if(labels & show_group)geom_text(aes(label = group), show.legend = FALSE)} +
      {if(labels & !show_group)geom_text(aes(label = sample), show.legend = FALSE)}
  }
  
  return(list(mds_values = mds, mds_plot = mds_plot))
  
}

