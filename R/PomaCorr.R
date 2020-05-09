
#' Correlation Analysis
#'
#' @description This function returns different correlation plots in both network and 2D formats and a table with all pairwise correlations in the data.
#' 
#' @param data A MSnSet object. First `pData` column must be the subject group/type.
#' @param shape Visualization shape of correlation matrix. Allowed values are "square" (default) and "circle".
#' @param type Character. Options are "full" (default), "lower" or "upper".
#' @param show_labels Logical indicating if correlation coefficient for each pair of features should be plotted or not (default = FALSE). Only recomended for a few number of features.
#' @param low Colour for low end of the gradient.
#' @param outline Colour for the outline of the gradient.
#' @param high Colour for high end of the gradient.
#' @param corr_type Type of graph to be made with correlation matrix. Options are "cor" (for global correlations), "pcor" (for partial correlations) and "glasso" (for gaussian graphical models). Default is "cor". See `glasso` R package for the last option.
#' @param threshold Numeric indicatin correlation coefficient. Edges with absolute weight below this value are will be removed from the network. 
#' @param rho Only if "corr_type" is set to "glasso". This parameter indicates the regularization parameter for lasso (rho = 0 means no regularization).
#' @param edge_labels Logical indicating if correlation coefficient will be plotted or not. Default is FALSE.
#' 
#' @export
#'
#' @return A list with the results.
#' @references Jerome Friedman, Trevor Hastie and Rob Tibshirani (2019). glasso: Graphical Lasso: Estimation of Gaussian Graphical Models. R package version 1.11. https://CRAN.R-project.org/package=glasso
#' @references Sacha Epskamp, Angelique O. J. Cramer, Lourens J. Waldorp, Verena D. Schmittmann, Denny Borsboom (2012). qgraph: Network Visualizations of Relationships in Psychometric Data. Journal of Statistical Software, 48(4), 1-18. URL http://www.jstatsoft.org/v48/i04/.
#' @author Pol Castellano-Escuder
#'
#' @importFrom ggplot2 theme_bw
#' @importFrom crayon red
#' @importFrom clisymbols symbol
#' @importFrom Biobase exprs
#' @importFrom ggcorrplot ggcorrplot
#' @importFrom qgraph qgraph
#' @importFrom glasso glasso
PomaCorr <- function(data,
                     shape = "square",
                     type = "full",
                     show_corr = FALSE,
                     low = "#336B87",
                     outline = "white",
                     high = "#EA8620",
                     corr_type = "cor",
                     threshold = 0.4,
                     rho = 0.7,
                     edge_labels = FALSE){
  
  if (missing(data)) {
    stop(crayon::red(clisymbols::symbol$cross, "data argument is empty!"))
  }
  if(!(class(data) == "MSnSet")){
    stop(paste0(crayon::red(clisymbols::symbol$cross, "data is not a MSnSet object."), 
                " \nSee POMA::PomaMSnSetClass or MSnbase::MSnSet"))
  }
  if (!(shape %in% c("square", "circle"))) {
    stop(crayon::red(clisymbols::symbol$cross, "Incorrect value for shape argument!"))
  }
  if (!(type %in% c("full", "lower"))) {
    stop(crayon::red(clisymbols::symbol$cross, "Incorrect value for type argument!"))
  }
  if (!(corr_type %in% c("cor", "pcor", "glasso"))) {
    stop(crayon::red(clisymbols::symbol$cross, "Incorrect value for corr_type argument!"))
  }
  
  total <- t(Biobase::exprs(data))
  
  ## All pairwise correlations
  cor_matrix <- cor(total)

  correlations <- cor_matrix
  correlations[lower.tri(correlations, diag = TRUE)] <- NA
  correlations <- as.data.frame(as.table(correlations))
  correlations <- na.omit(correlations)
  correlations <- correlations[with(correlations, order(-Freq)), ]
  colnames(correlations)[3] <- "corr"

  ## Corrplot
  my_cols <- c(low, outline, high)
  
  corrplot <- ggcorrplot(cor_matrix, method = shape, lab = show_corr, type = type,
                         ggtheme = ggplot2::theme_bw, colors = my_cols, legend.title = "Correlation")
  
  ## Networks
  node_names <- rownames(cor_matrix)
  
  if(corr_type != "glasso"){
    
    graph <- qgraph::qgraph(cor_matrix, graph = corr_type, layout = "spring", threshold = threshold, 
                            labels = substr(node_names, start = 1, stop = 5), edge.labels = edge_labels)
  } else{
    
    data_glasso <- glasso::glasso(cor_matrix, rho = rho)
    graph <- qgraph::qgraph(data_glasso, layout = "spring", threshold = threshold, labels = substr(node_names, start = 1, stop = 5),
                            edge.labels = edge_labels)
  }
  
  return(list(correlations = correlations, corrplot = corrplot, graph = graph))
  
}

