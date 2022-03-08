
#' Correlation Analysis
#'
#' @description This function returns different correlation plots and a table with all pairwise correlations in the data.
#' 
#' @param data A SummarizedExperiment object. First `colData` column must be the subject group/type.
#' @param method Character indicating which correlation coefficient has to be computed. Options are "pearson" (default), "kendall" and "spearman".
#' @param low Color for low end of the gradient in corrplot.
#' @param outline Color for the outline of the gradient in corrplot.
#' @param high Color for high end of the gradient in corrplot.
#' @param label_size Numeric indicating label size in corrplot.
#' @param corr_type Type of correlation network. Options are "cor" (for global correlations) and "glasso" (for gaussian graphical model). Default is "cor". See `glasso` R package for the second option.
#' @param coeff Numeric indicating correlation coefficient. Edges with absolute weight below this value will be removed from the network. If "corr_type" is set to "glasso", this parameter indicates the regularization parameter for lasso (rho = 0 means no regularization). See `glasso::glasso()`.
#' 
#' @export
#'
#' @return A list with the results.
#' @references Jerome Friedman, Trevor Hastie and Rob Tibshirani (2019). glasso: Graphical Lasso: Estimation of Gaussian Graphical Models. R package version 1.11. https://CRAN.R-project.org/package=glasso
#' @author Pol Castellano-Escuder
#'
#' @importFrom ggplot2 theme_bw
#' @importFrom dplyr filter rename as_tibble arrange desc mutate
#' @importFrom tidyr drop_na pivot_longer replace_na
#' @importFrom tibble rownames_to_column
#' @importFrom magrittr %>%
#' @importFrom SummarizedExperiment assay colData
#' @importFrom glasso glasso
#' 
#' @examples
#' data("st000284")
#' 
#' # Pearson correlation
#' PomaCorr(st000284)$correlations
#' 
#' ## Gaussian graphical model
#' # library(ggraph)
#' # PomaCorr(st000284, corr_type = "glasso")
PomaCorr <- function(data,
                     method = "pearson",
                     low = "red",
                     outline = "white",
                     high = "blue",
                     label_size = 12,
                     corr_type = "cor",
                     coeff = 0.7){
  
  if (missing(data)) {
    stop("data argument is empty!")
  }
  if(!is(data[1], "SummarizedExperiment")){
    stop("data is not a SummarizedExperiment object. \nSee POMA::PomaSummarizedExperiment or SummarizedExperiment::SummarizedExperiment")
  }
  if (!(corr_type %in% c("cor", "glasso"))) {
    stop("Incorrect value for corr_type argument!")
  }
  if (coeff > 1 | coeff < 0) {
    stop("coeff must be a number between 0 and 1...")
  }
  if (!(method %in% c("pearson", "kendall", "spearman"))) {
    stop("Incorrect value for method argument!")
  }
  
  e <- t(SummarizedExperiment::assay(data))
  
  cor_pmat <- function(x, ...) {

    mat <- as.matrix(x)
    n <- ncol(mat)
    p_mat <- matrix(NA, n, n)
    diag(p_mat) <- 0

    for (i in 1:(n - 1)) {
      for (j in (i + 1):n) {
        tmp <- stats::cor.test(mat[, i], mat[, j], method = method, ...)
        p_mat[i, j] <- p_mat[j, i] <- tmp$p.value
      }
    }
    
    colnames(p_mat) <- rownames(p_mat) <- colnames(mat)
    p_mat
  }
  
  flattenCorrMatrix <- function(cormat, pmat) {
    ut <- upper.tri(cormat)
    data.frame(
      row = rownames(cormat)[row(cormat)[ut]],
      column = rownames(cormat)[col(cormat)[ut]],
      cor  =(cormat)[ut],
      p = pmat[ut]
    )
  }
  
  cor_matrix <- cor(e, method = method)
  cor_pval <- cor_pmat(e)
  
  correlations <- flattenCorrMatrix(cor_matrix, cor_pval) %>% 
    dplyr::rename(feature1 = row, feature2 = column, R = cor, pvalue = p) %>% 
    tidyr::drop_na() %>% 	
    dplyr::mutate(FDR = p.adjust(pvalue, method = "fdr")) %>% 
    dplyr::arrange(desc(R)) %>% 
    dplyr::as_tibble()

  # corrplot

  cor_matrix_plot <- cor_matrix %>% 
    as.data.frame() %>% 
    tibble::rownames_to_column("feature1") %>% 
    tidyr::pivot_longer(cols = -feature1) %>% 
    dplyr::rename(feature2 = name, R = value) %>% 
    dplyr::mutate(R = tidyr::replace_na(R, 0)) %>% 
    dplyr::arrange(desc(R))
 
  my_cols <- c(low, outline, high)
  
  corrplot <- ggplot(cor_matrix_plot, aes_string(x = "feature1", y = "feature2", fill = "R")) +
    geom_tile(color = outline) +
    scale_fill_gradient2(
      low = my_cols[1],
      high = my_cols[3],
      mid = my_cols[2],
      midpoint = 0,
      limit = c(-1, 1),
      space = "Lab") +
    theme_bw() +
    theme(
      axis.text.x = element_text(
        angle = 45,
        vjust = 1,
        size = label_size,
        hjust = 1),
      axis.text.y = element_text(size = label_size),
      axis.title.x = element_blank(),
      axis.title.y = element_blank()
    ) +
    coord_fixed()
  
  # graph
  
  if(!(require("ggraph", character.only = TRUE))) {
    
    return(list(correlations = correlations, 
                corrplot = corrplot))
    
    warning("Package 'ggraph' is required for this function\nUse 'install.packages('ggraph')'")
  } else {
    
    if(corr_type != "glasso"){
      
      graph_table <- correlations %>% 
        filter(abs(R) >= coeff)
      
      if (nrow(graph_table) < 1) {
        stop("There are no feature pairs with selected coeff. Try with a lower value...")
      }
      
      graph <- ggraph(graph_table, layout = "fr") +
        geom_edge_link(aes(edge_alpha = abs(R), edge_width = abs(R), color = R)) +
        guides(edge_alpha = "none", edge_width = "none") +
        scale_edge_colour_gradientn(limits = c(-1, 1), colors = c("firebrick2", "dodgerblue2")) +
        geom_node_point(color = "white", size = 5) +
        geom_node_label(aes(label = name), repel = FALSE) +
        theme_graph()
      
      return(list(correlations = correlations, 
                  corrplot = corrplot, 
                  graph = graph))
      
    } else {
      
      data_glasso <- glasso::glasso(cor_matrix, rho = coeff)$w
      
      rownames(data_glasso) <- rownames(cor_matrix)
      colnames(data_glasso) <- colnames(cor_matrix)
      data_glasso[lower.tri(data_glasso, diag = TRUE)] <- NA
      data_glasso <- as.data.frame(as.table(data_glasso))
      data_glasso <- na.omit(data_glasso)
      data_glasso <- data_glasso[with(data_glasso, order(-Freq)), ]
      data_glasso <- data_glasso %>% 
        dplyr::rename(feature1 = 1, feature2 = 2, EstimatedCorr = 3) %>% 
        dplyr::as_tibble()
      
      graph_table <- data_glasso %>% 
        filter(EstimatedCorr != 0)
      
      if (nrow(graph_table) < 1) {
        stop("There aren't feature pairs with selected coeff. Try with a lower value...")
      }
      
      graph <- ggraph(graph_table, layout = "fr") +
        geom_edge_link(aes(edge_alpha = abs(EstimatedCorr), edge_width = abs(EstimatedCorr), color = EstimatedCorr)) +
        guides(edge_alpha = "none", edge_width = "none") +
        scale_edge_colour_gradientn(limits = c(-1, 1), colors = c("firebrick2", "dodgerblue2")) +
        geom_node_point(color = "white", size = 5) +
        geom_node_label(aes(label = name), repel = FALSE) +
        theme_graph()
      
      return(list(correlations = correlations, 
                  corrplot = corrplot, 
                  graph = graph, 
                  data_glasso = data_glasso))
      
    }
  }
}

