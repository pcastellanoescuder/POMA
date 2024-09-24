
#' Rank Product/Rank Sum Analysis
#'
#' @description `PomaRankProd` performs the Rank Product (or Rank Sum) method to identify differentially expressed genes.
#'
#' @param data A `SummarizedExperiment` object.
#' @param logged Logical. Indicates if data should be log transformed first.
#' @param paired Numeric. Indicates the number of random pairs generated in the function, if set to NA (default), the odd integer closer to the square of the number of replicates is used.
#' @param cutoff Numeric. Indicates the pfp/pvalue threshold value used to select features. Default is 1 to include all features.
#' @param method Character. Indicates the method to identify features. "pfp" uses percentage of false prediction, which is a default setting. "pval" uses p-values which is less stringent than pfp.
#'
#' @export
#'
#' @return A `list` with the results. Objects in the list are `up_regulated` (tibble) and `down_regulated` (tibble).
#' 
#' @references Breitling, R., Armengaud, P., Amtmann, A., and Herzyk, P.(2004) Rank Products: A simple, yet powerful, new method to detect differentially regulated genes in replicated microarray experiments, FEBS Letter, 57383-92
#' @references Hong, F., Breitling, R., McEntee, W.C., Wittner, B.S., Nemhauser, J.L., Chory, J. (2006). RankProd: a bioconductor package for detecting differentially expressed genes in meta-analysis Bioinformatics. 22(22):2825-2827
#' @references Del Carratore, F., Jankevics, A., Eisinga, R., Heskes, T., Hong, F. & Breitling, R. (2017). RankProd 2.0: a refactored Bioconductor package for detecting differentially expressed features in molecular profiling datasets. Bioinformatics. 33(17):2774-2775
#'
#' @author Pol Castellano-Escuder
#'
#' @importFrom magrittr %>%
#' 
#' @examples 
#' data <- POMA::st000336 %>% # Example SummarizedExperiment object included in POMA
#'   PomaImpute()
#' 
#' ## Output is a list with objects `up_regulated` (tibble with up regulated features) and `down_regulated` (tibble with down regulated features) 
#' ## Perform on no-scaled object to avoid negative values
#' data %>% 
#'   PomaRankProd(method = "pfp")
PomaRankProd <- function(data,
                         logged = TRUE,
                         paired = NA,
                         cutoff = 1,
                         method = "pfp") {

  if(!is(data, "SummarizedExperiment")){
    stop("data is not a SummarizedExperiment object. \nSee POMA::PomaCreateObject or SummarizedExperiment::SummarizedExperiment")
  }
  if (!(method %in% c("pfp", "pval"))) {
    stop("Incorrect value for method argument")
  }
  if (sum(apply(t(SummarizedExperiment::assay(data)), 2, function(x){sum(x < 0, na.rm = TRUE)})) != 0) {
    stop("Negative values nor allowed")
  }
  
  group_factor <- SummarizedExperiment::colData(data)[,1]
  
  if (!is.factor(group_factor)) {
    stop("Grouping factor must be a factor (first column of the metadata file)")
  }
  if (length(table(group_factor)[table(group_factor) != 0]) != 2) {
    stop("Grouping factor must have exactly 2 levels (first column of the metadata file)")
  }
  
  data_class <- as.numeric(ifelse(group_factor == levels(group_factor)[1], 0, 1))
  
  class1 <- levels(as.factor(group_factor))[1]
  class2 <- levels(as.factor(group_factor))[2]

  capture.output({
    RP <- RankProd::RankProducts(SummarizedExperiment::assay(data), 
                                 data_class, 
                                 logged = logged, 
                                 na.rm = TRUE, 
                                 plot = FALSE,
                                 RandomPairs = paired,
                                 rand = 123,
                                 gene.names = rownames(data))
    
    top_rank <- RankProd::topGene(RP, 
                                  cutoff = cutoff, 
                                  method = method,
                                  logged = logged, 
                                  logbase = 2,
                                  gene.names = rownames(data))
  }, file = "/dev/null")
          
  one <- as.data.frame(top_rank$Table1)
  two <- as.data.frame(top_rank$Table2)
  
  if (nrow(one) != 0){
    
    one <- one %>% 
      tibble::rownames_to_column("feature") %>% 
      dplyr::rename(rp_rsum = 3,
                    pvalue = P.value,
                    feature_index = gene.index) %>% 
      dplyr::as_tibble()
    
    colnames(one)[4] <- paste0("FC_", class1, "_", class2)
  }
  
  if (nrow(two) != 0) {
    
    two <- two %>% 
      tibble::rownames_to_column("feature") %>% 
      dplyr::rename(rp_rsum = 3,
                    pvalue = P.value,
                    feature_index = gene.index) %>% 
      dplyr::as_tibble()
    
    colnames(two)[4] <- paste0("FC_", class1, "_", class2)
  }

  # Plot
  # pfp <- as.matrix(RP$pfp)
  # 
  # if (is.null(RP$RPs)) {
  #   RP1 <- as.matrix(RP$RSs)
  #   rank <- as.matrix(RP$RSrank)
  # }
  # 
  # if (!is.null(RP$RPs)){
  #   RP1 <- as.matrix(RP$RPs)
  #   rank <- as.matrix(RP$RPrank)
  # }
  # 
  # ind1 <- which(!is.na(RP1[, 1]))
  # ind2 <- which(!is.na(RP1[, 2]))
  # ind3 <- append(ind1, ind2)
  # ind3 <- unique(ind3)
  # RP.sort.upin2 <- sort(RP1[ind1, 1], index.return = TRUE)
  # RP.sort.downin2 <- sort(RP1[ind2, 2], index.return = TRUE)
  # pfp1 <- pfp[ind1, 1]
  # pfp2 <- pfp[ind2, 2]
  # rank1 <- rank[ind1, 1]
  # rank2 <- rank[ind2, 2]
  # 
  # rp_plot <- data.frame(rank1 = rank1, rank2 = rank2, pfp1 = pfp1 ,  pfp2 = pfp2)
  # 
  # plot1 <- ggplot2::ggplot(rp_plot, ggplot2::aes(x = rank1, y = pfp1)) +
  #   ggplot2::geom_point(size = 1.5, alpha=0.9) +
  #   theme_poma() +
  #   ggplot2::labs(x = "Number of identified features",
  #                 y = "Estimated PFP",
  #                 title = paste0("Up-regulated features in ", class2))
  # 
  # plot2 <- ggplot2::ggplot(rp_plot, ggplot2::aes(x = rank2, y = pfp2)) +
  #   ggplot2::geom_point(size = 1.5, alpha=0.9) +
  #   theme_poma() +
  #   ggplot2::labs(x = "Number of identified features",
  #                 y = "Estimated PFP",
  #                 title = paste0("Down-regulated features in ", class2))

  return(list(up_regulated = one,
              down_regulated = two
              # up_regulated_plot = plot1,
              # down_regulated_plot = plot2
              ))
}

