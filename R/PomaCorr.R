
#' Correlation P-Values
#'
#' Compute correlation p-values.
#' 
#' @param x A data matrix.
#' @param method Character indicating which correlation coefficient has to be computed. Options are "pearson" (default), "kendall" and "spearman".
cor_pmat <- function(x, method) {
  
  mat <- as.matrix(x)
  n <- ncol(mat)
  p_mat <- matrix(NA, n, n)
  diag(p_mat) <- 0
  
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- stats::cor.test(mat[, i], mat[, j], method = method)
      p_mat[i, j] <- p_mat[j, i] <- tmp$p.value
    }
  }
  
  colnames(p_mat) <- rownames(p_mat) <- colnames(mat)
  p_mat
}

#' Flatten Correlation Matrix
#' 
#' @param cormat Output from `cor`.
#' @param pmat Output from `cor_pmat`.
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor = (cormat)[ut],
    p = pmat[ut]
  )
}

#' Correlation Analysis
#'
#' @description `PomaCorr` computes all pairwise correlations in the data.
#' 
#' @param data A `SummarizedExperiment` object.
#' @param method Character. Indicates which correlation coefficient has to be computed. Options are "pearson" (default), "kendall" and "spearman".
#' @param label_size Numeric. Indicates plot label size.
#' @param theme_params List. Indicating `theme_poma` parameters.
#' 
#' @export
#'
#' @return A `list` with the results.
#' @references Jerome Friedman, Trevor Hastie and Rob Tibshirani (2019). glasso: Graphical Lasso: Estimation of Gaussian Graphical Models. R package version 1.11. https://CRAN.R-project.org/package=glasso
#' @author Pol Castellano-Escuder
#' 
#' @importFrom magrittr %>%
#' 
#' @examples
#' ## Output is a list with objects `correlations` (tibble) and `corrplot` (ggplot2 object)
#' data <- POMA::st000284 # Example SummarizedExperiment object included in POMA
#' 
#' data %>% 
#'   PomaCorr(method = "pearson")
PomaCorr <- function(data,
                     method = "pearson",
                     label_size = 8,
                     theme_params = list()) {
  
  if (missing(data)) {
    stop("data argument is empty!")
  }
  if (!is(data, "SummarizedExperiment")){
    stop("data is not a SummarizedExperiment object. \nSee POMA::PomaCreateObject or SummarizedExperiment::SummarizedExperiment")
  }
  if (!(method %in% c("pearson", "kendall", "spearman"))) {
    stop("Incorrect value for method argument")
  }
  
  to_correlation <- t(SummarizedExperiment::assay(data))
  
  cor_matrix <- cor(to_correlation, method = method)
  cor_pval <- cor_pmat(to_correlation, method = method)
  
  correlations <- flattenCorrMatrix(cor_matrix, cor_pval) %>% 
    dplyr::rename(feature1 = row, feature2 = column, corr = cor, pvalue = p) %>% 
    tidyr::drop_na() %>% 	
    dplyr::mutate(FDR = p.adjust(pvalue, method = "fdr")) %>% 
    dplyr::arrange(dplyr::desc(corr)) %>% 
    dplyr::as_tibble()

  # corrplot
  cor_matrix_plot <- cor_matrix %>% 
    as.data.frame() %>% 
    tibble::rownames_to_column("feature1") %>% 
    tidyr::pivot_longer(cols = -feature1) %>% 
    dplyr::rename(feature2 = name, corr = value) %>% 
    dplyr::mutate(corr = tidyr::replace_na(corr, 0)) %>% 
    dplyr::arrange(dplyr::desc(corr))
  
  corrplot <- ggplot2::ggplot(cor_matrix_plot, ggplot2::aes_string(x = "feature1", y = "feature2", fill = "corr")) +
    ggplot2::geom_tile() +
    theme_poma(axis_x_rotate = TRUE, axistitle = "none", base_size = label_size) +
    ggplot2::scale_fill_gradient2(
      low = "red",
      high = "blue",
      mid = "white",
      midpoint = 0,
      limit = c(-1, 1),
      space = "Lab")
  
  # graph
  # if(!(requireNamespace("ggraph", character.only = TRUE))) {
  #   
  #   return(list(correlations = correlations, 
  #               corrplot = corrplot))
  #   
  #   warning("Package 'ggraph' is required for this function\nUse 'install.packages('ggraph')'")
  # } else {
  #   
  #   if(corr_type != "glasso"){
  #     
  #     graph_table <- correlations %>% 
  #       dplyr::filter(abs(corr) >= coeff)
  #     
  #     if (nrow(graph_table) < 1) {
  #       stop("There are no feature pairs with selected coeff. Use a lower value...")
  #     }
  #     
  #     graph <- ggraph::ggraph(graph_table, layout = "fr") +
  #       ggraph::geom_edge_link(ggplot2::aes(edge_alpha = abs(corr), edge_width = abs(corr), color = corr)) +
  #       ggplot2::guides(edge_alpha = "none", edge_width = "none") +
  #       ggraph::scale_edge_colour_gradientn(limits = c(-1, 1), colors = c("firebrick2", "dodgerblue2")) +
  #       ggraph::geom_node_point(color = "white", size = 2) +
  #       ggraph::geom_node_text(ggplot2::aes(label = name), repel = FALSE) +
  #       ggraph::theme_graph()
  #     
  #     return(list(correlations = correlations, 
  #                 corrplot = corrplot, 
  #                 graph = graph))
  #     
  #   } else {
  #     
  #     data_glasso <- glasso::glasso(cor_matrix, rho = coeff)$w
  #     
  #     rownames(data_glasso) <- rownames(cor_matrix)
  #     colnames(data_glasso) <- colnames(cor_matrix)
  #     data_glasso[lower.tri(data_glasso, diag = TRUE)] <- NA
  #     data_glasso <- as.data.frame(as.table(data_glasso))
  #     data_glasso <- na.omit(data_glasso)
  #     data_glasso <- data_glasso[with(data_glasso, order(-Freq)), ]
  #     data_glasso <- data_glasso %>% 
  #       dplyr::rename(feature1 = 1, feature2 = 2, est_corr = 3) %>% 
  #       dplyr::as_tibble()
  #     
  #     graph_table <- data_glasso %>% 
  #       dplyr::filter(est_corr != 0)
  #     
  #     if (nrow(graph_table) < 1) {
  #       stop("There aren't feature pairs with selected coeff. Try with a lower value...")
  #     }
  #     
  #     graph <- ggraph::ggraph(graph_table, layout = "fr") +
  #       ggraph::geom_edge_link(ggplot2::aes(edge_alpha = abs(est_corr), edge_width = abs(est_corr), color = est_corr)) +
  #       ggplot2::guides(edge_alpha = "none", edge_width = "none") +
  #       ggraph::scale_edge_colour_gradientn(limits = c(-1, 1), colors = c("firebrick2", "dodgerblue2")) +
  #       ggraph::geom_node_point(color = "white", size = 2) +
  #       ggraph::geom_node_text(ggplot2::aes(label = name), repel = FALSE) +
  #       ggraph::theme_graph()
  #     
  #     return(list(correlations = correlations, 
  #                 corrplot = corrplot, 
  #                 graph = graph, 
  #                 data_glasso = data_glasso))
  #     
  #   }
  # }
  
  return(list(correlations = correlations,
                  corrplot = corrplot)
         )
}

