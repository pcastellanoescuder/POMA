
#' Correlation Analysis
#'
#' @description This function returns different correlation plots (correlogram and network plots) and a table with all pairwise correlations in the data.
#' 
#' @param data A MSnSet object. First `pData` column must be the subject group/type.
#' @param method Character indicating which correlation coefficient has to be computed. Options are "pearson" (default), "kendall" and "spearman".
#' @param shape Character ingicating shape of correlogram. Options are "square" (default) and "circle".
#' @param type Character indicating type of correlogram. Options are "full" (default), "lower" or "upper".
#' @param show_corr Logical indicating if correlation coefficient for each pair of features should be plotted in correlogram or not (default = FALSE). Only recomended for a low number of features.
#' @param low Colour for low end of the gradient in correlogram.
#' @param outline Colour for the outline of the gradient in correlogram.
#' @param high Colour for high end of the gradient in correlogram.
#' @param label_size Numeric indicating label size in correlogram.
#' @param corr_type Type of network to be made with correlation matrix. Options are "cor" (for global correlations) and "glasso" (for gaussian graphical model). Default is "cor". See `glasso` R package for the second option.
#' @param coeff Numeric indicatin correlation coefficient. Edges with absolute weight below this value will be removed from the network. If "corr_type" is set to "glasso", this parameter indicates the regularization parameter for lasso (rho = 0 means no regularization). See `glasso::glasso()`.
#' 
#' @export
#'
#' @return A list with the results.
#' @references Jerome Friedman, Trevor Hastie and Rob Tibshirani (2019). glasso: Graphical Lasso: Estimation of Gaussian Graphical Models. R package version 1.11. https://CRAN.R-project.org/package=glasso
#' @author Pol Castellano-Escuder
#'
#' @import ggraph
#' @importFrom ggplot2 theme_bw
#' @importFrom dplyr filter
#' @importFrom magrittr %>%
#' @importFrom crayon red
#' @importFrom clisymbols symbol
#' @importFrom MSnbase exprs
#' @importFrom ggcorrplot ggcorrplot
#' @importFrom glasso glasso
#' 
#' @examples 
#' data("st000284")
#' 
#' # pearson correlation
#' PomaCorr(st000284)$correlations
#' PomaCorr(st000284)$corrplot
#' 
#' # gaussian graphical model
#' # library(ggraph)
#' # PomaCorr(st000284, corr_type = "glasso")
PomaCorr <- function(data,
                     method = "pearson",
                     shape = "square",
                     type = "full",
                     show_corr = FALSE,
                     low = "#336B87",
                     outline = "white",
                     high = "#EA8620",
                     label_size = 12,
                     corr_type = "cor",
                     coeff = 0.7){
  
  if (missing(data)) {
    stop(crayon::red(clisymbols::symbol$cross, "data argument is empty!"))
  }
  if(!is(data[1], "MSnSet")){
    stop(paste0(crayon::red(clisymbols::symbol$cross, "data is not a MSnSet object."), 
                " \nSee POMA::PomaMSnSetClass or MSnbase::MSnSet"))
  }
  if (!(shape %in% c("square", "circle"))) {
    stop(crayon::red(clisymbols::symbol$cross, "Incorrect value for shape argument!"))
  }
  if (!(type %in% c("full", "lower", "upper"))) {
    stop(crayon::red(clisymbols::symbol$cross, "Incorrect value for type argument!"))
  }
  if (!(corr_type %in% c("cor", "glasso"))) {
    stop(crayon::red(clisymbols::symbol$cross, "Incorrect value for corr_type argument!"))
  }
  if (coeff > 1 | coeff < 0) {
    stop(crayon::red(clisymbols::symbol$cross, "coeff must be a number between 0 and 1..."))
  }
  if (!(method %in% c("pearson", "kendall", "spearman"))) {
    stop(crayon::red(clisymbols::symbol$cross, "Incorrect value for method argument!"))
  }
  
  total <- t(MSnbase::exprs(data))
  cor_matrix <- cor(total, method = method)
  
  ## Pairwise correlations
  correlations <- cor_matrix
  correlations[lower.tri(correlations, diag = TRUE)] <- NA
  correlations <- as.data.frame(as.table(correlations))
  correlations <- na.omit(correlations)
  correlations <- correlations[with(correlations, order(-Freq)), ]
  colnames(correlations)[3] <- "corr"

  ## Corrplot
  my_cols <- c(low, outline, high)
  
  corrplot <- ggcorrplot(cor_matrix, method = shape, lab = show_corr, type = type,
                         ggtheme = ggplot2::theme_bw, colors = my_cols, 
                         legend.title = "Correlation", tl.cex = label_size)
  
  ## Networks
  if(corr_type != "glasso"){
    
    graph_table <- correlations %>% 
      filter(abs(corr) >= coeff)
    
    if (nrow(graph_table) < 1) {
      stop(crayon::red(clisymbols::symbol$cross, "There are no feature pairs with selected coeff. Try with a lower value..."))
    }
    
    graph <- ggraph(graph_table, layout = "fr") +
      geom_edge_link(aes(edge_alpha = abs(corr), edge_width = abs(corr), color = corr)) +
      guides(edge_alpha = "none", edge_width = "none") +
      scale_edge_colour_gradientn(limits = c(-1, 1), colors = c("firebrick2", "dodgerblue2")) +
      geom_node_point(color = "white", size = 5) +
      geom_node_label(aes(label = name), repel = FALSE) +
      theme_graph()
    
  } else {
    
    data_glasso <- glasso::glasso(cor_matrix, rho = coeff)$w
    
    rownames(data_glasso) <- rownames(cor_matrix)
    colnames(data_glasso) <- colnames(cor_matrix)
    data_glasso[lower.tri(data_glasso, diag = TRUE)] <- NA
    data_glasso <- as.data.frame(as.table(data_glasso))
    data_glasso <- na.omit(data_glasso)
    data_glasso <- data_glasso[with(data_glasso, order(-Freq)), ]
    colnames(data_glasso)[3] <- "EstimatedCorr"
    
    graph_table <- data_glasso %>% 
      filter(EstimatedCorr != 0)
    
    if (nrow(graph_table) < 1) {
      stop(crayon::red(clisymbols::symbol$cross, "There are no feature pairs with selected coeff. Try with a lower value..."))
    }
    
    graph <- ggraph(graph_table, layout = "fr") +
      geom_edge_link(aes(edge_alpha = abs(EstimatedCorr), edge_width = abs(EstimatedCorr), color = EstimatedCorr)) +
      guides(edge_alpha = "none", edge_width = "none") +
      scale_edge_colour_gradientn(limits = c(-1, 1), colors = c("firebrick2", "dodgerblue2")) +
      geom_node_point(color = "white", size = 5) +
      geom_node_label(aes(label = name), repel = FALSE) +
      theme_graph()

  }
  
  if(corr_type != "glasso"){
    return(list(correlations = correlations, corrplot = corrplot, graph = graph))
  } else{
    return(list(correlations = correlations, corrplot = corrplot, graph = graph, data_glasso = data_glasso))
  }
  
}

