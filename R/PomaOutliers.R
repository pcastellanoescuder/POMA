
#' Remove and Analyze Outliers
#'
#' @description blabla
#'
#' @param data A MSnSet object. First `pData` column must be the subject group/type.
#' @param do Action to do. Options are "clean" (to remove detected outliers) and "analyze" (to analyze data outliers). Note that the output of this function will be different depending on this parameter.
#' @param method blabla
#' @param type blabla
#' @param iqr blabla
#' @param labels blabla
#'
#' @export
#'
#' @return A MSnSet object with cleaned data or different plots for the detailed analysis of outliers (depending on "do" parameter).
#' @author Pol Castellano-Escuder
#'
#' @import ggplot2
#' @importFrom dplyr select filter do group_by mutate rename
#' @importFrom tibble rownames_to_column
#' @importFrom magrittr %>%
#' @importFrom ggrepel geom_label_repel geom_text_repel
#' @importFrom crayon red
#' @importFrom clisymbols symbol
#' @importFrom Biobase pData exprs sampleNames varLabels
#' @importFrom vegan betadisper
PomaOutliers <- function(data, # nocov start
                         do = "clean",
                         method = "euclidean",
                         type = "median",
                         iqr = 1.5,
                         labels = FALSE){
  
  Biobase::varLabels(data)[1] <- "Group"
  groups <- Biobase::pData(data)$Group
  names <- Biobase::sampleNames(data)
  to_outliers <- t(Biobase::exprs(data))
  
  ##
  
  dd <- dist(to_outliers, method = method)
  distances <- vegan::betadisper(dd, groups, type = type, bias.adjust = FALSE, sqrt.dist = FALSE, add = FALSE)
  detect_outliers <- data.frame(distances = distances$distances, Groups = distances$group) 
  
  limit <- data.frame(aggregate(detect_outliers$distances, list(detect_outliers$Groups), function(x) {quantile(x, 0.75) + iqr * IQR(x)}))
  colnames(limit)[1] <- "Groups"
  
  detect_outliers <- merge(detect_outliers, limit, by = "Groups")
  detect_outliers <- detect_outliers %>% 
    mutate(outlier = as.factor(ifelse(distances > x, 1, 0)),
           sample = names)
  
  final_outliers <- detect_outliers %>%
    filter(outlier == 1) %>%
    select(sample, Groups, distances, x) %>%
    rename(group = Groups,
           distance_to_centroid = distances,
           limit_distance = x)
  
  ##
  
  if(do == "analyze"){
    
    vectors <- data.frame(distances$vectors)
    centroids <- data.frame(distances$centroids)
    
    total_outliers <- vectors %>%
      mutate(Group = groups) 
    
    find_hull <- function(x){x[chull(x$PCoA1, x$PCoA2) ,]}
    
    hulls <- total_outliers %>% 
      group_by(Group) %>% 
      dplyr::do(find_hull(.))
    
    polygon_plot <- ggplot(total_outliers, aes(x = PCoA1, y = PCoA2)) +
      geom_polygon(data = hulls, alpha = 0.5, aes(fill = Group)) +
      geom_point(aes(shape = Group), size = 3) +
      geom_label(data = centroids, aes(x = PCoA1, y = PCoA2, label = rownames(centroids)), show.legend = F) +
      {if(labels)ggrepel::geom_text_repel(aes(label = rownames(total_outliers)))} +
      theme_bw()
    
    distance_boxplot <- ggplot(detect_outliers, aes(Groups, distances, fill = Groups)) +
      geom_boxplot() +
      ylab("Distance to group centroid") + 
      xlab("") +
      {if(labels)ggrepel::geom_label_repel(data = detect_outliers[detect_outliers$outlier == 1,], aes(label = sample), na.rm = TRUE, size = 4, show.legend = F)} +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    return(list(polygon_plot = polygon_plot, distance_boxplot = distance_boxplot, outliers = final_outliers))
    
  } else {
    
    target <- pData(data) %>% 
      rownames_to_column("sample") %>% 
      as.data.frame() %>%
      filter(!(sample %in% final_outliers$sample))

    e <- to_outliers %>%
      as.data.frame() %>%
      filter(!(rownames(.) %in% final_outliers$sample))
    
    ##
    
    dataCleaned <- PomaMSnSetClass(features = e, target = target)
    
    dataCleaned@processingData@processing <-
      c(data@processingData@processing,
        paste("Outliers removed (", method , " and ", type, "): ", date(), sep = ""))
    dataCleaned@processingData@cleaned <- TRUE
    dataCleaned@experimentData <- data@experimentData
    dataCleaned@qual <- data@qual
    
    if (validObject(dataCleaned))
      return(dataCleaned)
    
  }
  
} # nocov end

