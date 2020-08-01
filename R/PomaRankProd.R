
#' Rank Product/Rank Sum Analysis for Mass Spectrometry Data
#'
#' @description PomaRankProd() performs the Rank Product method to identify differential feature concentration/intensity.
#'
#' @param data A MSnSet object. First `pData` column must be the subject group/type.
#' @param logged If "TRUE" (default) data have been previously log transformed.
#' @param logbase Numerical. Base for log transformation.
#' @param paired Number of random pairs generated in the function, if set to NA (default), the odd integer closer to the square of the number of replicates is used.
#' @param cutoff The pfp/pvalue threshold value used to select features.
#' @param method If cutoff is provided, the method needs to be selected to identify features. "pfp" uses percentage of false prediction, which is a default setting. "pval" uses p-values which is less stringent than pfp.
#'
#' @export
#'
#' @return A list with all results for Rank Product analysis including tables and plots.
#' @references Breitling, R., Armengaud, P., Amtmann, A., and Herzyk, P.(2004) Rank Products: A simple, yet powerful, new method to detect differentially regulated genes in replicated microarray experiments, FEBS Letter, 57383-92
#' @author Pol Castellano-Escuder
#'
#' @importFrom RankProd RankProducts topGene
#' @import ggplot2
#' @importFrom crayon red blue
#' @importFrom clisymbols symbol
#' @importFrom Biobase varLabels pData exprs
PomaRankProd <- function(data,
                         logged = TRUE,
                         logbase = 2,
                         paired = NA,
                         cutoff = 0.05,
                         method = "pfp"){

  if (missing(data)) {
    stop(crayon::red(clisymbols::symbol$cross, "data argument is empty!"))
  }
  if(!is(data[1], "MSnSet")){
    stop(paste0(crayon::red(clisymbols::symbol$cross, "data is not a MSnSet object."), 
                " \nSee POMA::PomaMSnSetClass or MSnbase::MSnSet"))
  }
  if (!(method %in% c("pfp", "pval"))) {
    stop(crayon::red(clisymbols::symbol$cross, "Incorrect value for method argument!"))
  }
  if (missing(method)) {
    warning("method argument is empty! pfp method will be used")
  }

  Biobase::varLabels(data)[1] <- "Group"
  Group <- as.factor(Biobase::pData(data)$Group)

  if (length(levels(Group)) != 2) {
    stop(crayon::red(clisymbols::symbol$cross, "Data must have two groups..."))
  }

  data_class <- as.numeric(ifelse(Group == levels(Group)[1], 0, 1))
  
  class1 <- levels(as.factor(Group))[1]
  class2 <- levels(as.factor(Group))[2]

  RP <- RankProducts(Biobase::exprs(data), data_class, logged = logged, na.rm = TRUE, plot = FALSE,
                     RandomPairs = paired,
                     rand = 123,
                     gene.names = rownames(data))

  top_rank <- topGene(RP, cutoff = cutoff, method = method,
                      logged = logged, logbase = logbase,
                      gene.names = rownames(data))

  one <- as.data.frame(top_rank$Table1)
  two <- as.data.frame(top_rank$Table2)

  if(nrow(one) == 0 & nrow(two) == 0){
    stop(crayon::blue(clisymbols::symbol$info, "No significant features found..."))
  }
  
  colnames(one)[3] <- paste0("FC: ", class1, "/", class2)
  colnames(two)[3] <- paste0("FC: ", class1, "/", class2)

  one <- one %>% dplyr::select(-gene.index)
  two <- two %>% dplyr::select(-gene.index)

  #### PLOT

  pfp <- as.matrix(RP$pfp)

  ####

  if (is.null(RP$RPs)) {
    RP1 <- as.matrix(RP$RSs)
    rank <- as.matrix(RP$RSrank)
  }

  if (!is.null(RP$RPs)){
    RP1 <- as.matrix(RP$RPs)
    rank <- as.matrix(RP$RPrank)
  }

  ind1 <- which(!is.na(RP1[, 1]))
  ind2 <- which(!is.na(RP1[, 2]))
  ind3 <- append(ind1, ind2)
  ind3 <- unique(ind3)
  RP.sort.upin2 <- sort(RP1[ind1, 1], index.return = TRUE)
  RP.sort.downin2 <- sort(RP1[ind2, 2], index.return = TRUE)
  pfp1 <- pfp[ind1, 1]
  pfp2 <- pfp[ind2, 2]
  rank1 <- rank[ind1, 1]
  rank2 <- rank[ind2, 2]

  rp_plot <- data.frame(rank1 = rank1, rank2 = rank2, pfp1 = pfp1 ,  pfp2 = pfp2)

  plot1 <- ggplot(rp_plot, aes(x = rank1, y = pfp1)) +
    geom_point(size = 1.5, alpha=0.8) +
    theme_bw() +
    xlab("Number of identified features") +
    ylab("Estimated PFP") +
    ggtitle(paste0("Identification of Up-regulated features under class ", class2))

  plot2 <- ggplot(rp_plot, aes(x = rank2, y = pfp2)) +
    geom_point(size = 1.5, alpha=0.8) +
    theme_bw() +
    xlab("Number of identified features") +
    ylab("Estimated PFP") +
    ggtitle(paste0("Identification of Down-regulated features under class ", class2))

  return(list(upregulated = one,
              downregulated = two,
              Upregulated_RP_plot = plot1,
              Downregulated_RP_plot = plot2))

}

