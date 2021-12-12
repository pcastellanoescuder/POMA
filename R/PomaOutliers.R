
#' Remove and Analyze Outliers
#'
#' @description This function allows users to analyze outliers by different plots and remove them from an SummarizedExperiment object.
#'
#' @param data A SummarizedExperiment object. First `colData` column must be the subject group/type.
#' @param do Action to do. Options are "clean" (to remove detected outliers) and "analyze" (to analyze data outliers). Note that the output of this function will be different depending on this parameter.
#' @param method Distance measure method to perform MDS. Options are "euclidean", "maximum", "manhattan", "canberra" and "minkowski". See `?dist()`.
#' @param type Type of outliers analysis to perform. Options are "median" (default) and "centroid". See `vegan::betadisper`.
#' @param coef This value corresponds to the classical 1.5 in \eqn{Q3 + 1.5*IQR} formula to detect outliers. By changing this value, the permissiveness in outlier detection will change.
#' @param labels Logical indicating if sample IDs should to be plotted or not.
#'
#' @export
#'
#' @return A SummarizedExperiment object with cleaned data or different exploratory plots for the detailed analysis of outliers (depending on "do" parameter).
#' @author Pol Castellano-Escuder
#'
#' @import ggplot2
#' @importFrom dplyr select filter do group_by mutate rename
#' @importFrom tibble rownames_to_column
#' @importFrom magrittr %>%
#' @importFrom ggrepel geom_label_repel
#' @importFrom SummarizedExperiment assay colData
#' @importFrom vegan betadisper
#' 
#' @examples 
#' data("st000336")
#' 
#' # clean outliers
#' st000336 %>% 
#'   PomaImpute() %>%
#'   PomaNorm() %>%
#'   PomaOutliers()
#' 
#' # analyze outliers
#' st000336 %>% 
#'   PomaImpute() %>%
#'   PomaNorm() %>%
#'   PomaOutliers(do = "analyze")
PomaOutliers <- function(data,
                         do = "clean",
                         method = "euclidean",
                         type = "median",
                         coef = 1.5,
                         labels = FALSE){
  
  if (missing(data)) {
    stop("data argument is empty!")
  }
  if(!is(data[1], "SummarizedExperiment")){
    stop("data is not a SummarizedExperiment object. \nSee POMA::PomaSummarizedExperiment or SummarizedExperiment::SummarizedExperiment")
  }
  if (!(method %in% c("euclidean", "maximum", "manhattan", "canberra", "minkowski"))) {
    stop("Incorrect value for method argument!")
  }
  if (!(type %in% c("median", "centroid"))) {
    stop("Incorrect value for type argument!")
  }
  if (!(do %in% c("clean", "analyze"))) {
    stop("Incorrect value for do argument!")
  }
  
  groups <- SummarizedExperiment::colData(data)[,1]
  names <- rownames(SummarizedExperiment::colData(data))
  to_outliers <- t(SummarizedExperiment::assay(data))
  
  ##
  
  dd <- dist(to_outliers, method = method)
  distances <- vegan::betadisper(dd, groups, type = type, bias.adjust = FALSE, sqrt.dist = FALSE, add = FALSE)
  detect_outliers <- data.frame(distances = distances$distances, Groups = distances$group, sample = names) 
  
  limit <- data.frame(aggregate(detect_outliers$distances, list(detect_outliers$Groups), function(x) {quantile(x, 0.75) + coef * IQR(x)}))
  colnames(limit)[1] <- "Groups"
  
  detect_outliers <- merge(detect_outliers, limit, by = "Groups")
  detect_outliers <- detect_outliers %>% 
    mutate(out = as.factor(ifelse(distances > x, 1, 0)))
  
  final_outliers <- detect_outliers %>%
    filter(out == 1) %>%
    select(sample, Groups, distances, x) %>%
    rename(group = Groups,
           distance_to_centroid = distances,
           limit_distance = x)
  
  ##
  
  if(do == "analyze"){
    
    vectors <- data.frame(distances$vectors)
    centroids <- data.frame(distances$centroids)
    
    total_outliers <- vectors %>%
      rownames_to_column("sample") %>%
      mutate(Group = groups)
    
    find_hull <- function(x){x[chull(x$PCoA1, x$PCoA2) ,]}
    
    hulls <- total_outliers %>% 
      select(-sample) %>%
      group_by(Group) %>% 
      dplyr::do(find_hull(.))
    
    polygon_plot <- ggplot(total_outliers, aes(x = PCoA1, y = PCoA2)) +
      geom_polygon(data = hulls, alpha = 0.5, aes(fill = Group)) +
      {if(!labels)geom_point(aes(shape = Group), size = 3, alpha = 0.7)} +
      geom_label(data = centroids, aes(x = PCoA1, y = PCoA2, color = rownames(centroids), label = rownames(centroids)), show.legend = FALSE) +
      {if(labels)geom_text(aes(label = sample))} +
      theme_bw() +
      scale_fill_viridis_d()
    
    distance_boxplot <- ggplot(detect_outliers, aes(Groups, distances, fill = Groups)) +
      geom_boxplot(coef = coef, alpha = 0.8) +
      ylab("Distance to group centroid") + 
      xlab("") +
      {if(labels)ggrepel::geom_label_repel(data = detect_outliers[detect_outliers$out == 1,], aes(label = sample), na.rm = TRUE, size = 4, show.legend = FALSE)} +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      scale_fill_viridis_d()
    
    return(list(polygon_plot = polygon_plot, distance_boxplot = distance_boxplot, outliers = final_outliers))
    
  } else {
    
    target <- SummarizedExperiment::colData(data) %>% 
      as.data.frame() %>%
      rownames_to_column("sample") %>% 
      filter(!(sample %in% final_outliers$sample))

    e <- to_outliers %>%
      as.data.frame() %>%
      filter(!(rownames(.) %in% final_outliers$sample))
    
    ##
    
    dataCleaned <- PomaSummarizedExperiment(features = e, target = target)
    
    if (validObject(dataCleaned))
      return(dataCleaned)
    
  }
  
}

