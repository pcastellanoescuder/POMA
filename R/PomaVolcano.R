
#' Volcano Plot
#'
#' @description `PomaVolcano` creates a volcano plot from a given dataset. This function is designed to visualize the statistical significance (p-value) against the magnitude of change (log2 fold change) for features.
#'
#' @param data A data frame with at least three columns: feature names, effect size (e.g., logFC), and statistical significance values. These should be the first three columns of the data, in this order.
#' @param pval_cutoff Numeric. Specifies the p-value threshold for significance in the volcano plot. The default is set to 0.05. This parameter determines the horizontal line in the plot indicating the threshold for statistical significance.
#' @param log2fc_cutoff Numeric. Specifies the log2 fold change cutoff for the volcano plot. If not provided, the cutoff is set to the 75th percentile of the absolute log2 fold changes in the data. This parameter determines the vertical lines in the plot indicating the magnitude of change threshold.
#' @param labels Logical. Indicates whether to plot labels for significant features. 
#' @param x_label Character. Custom label for the x-axis.
#' @param y_label Character. Custom label for the y-axis.
#'
#' @export
#' 
#' @return A `ggplot` object representing the volcano plot.
#' @author Pol Castellano-Escuder
#' 
#' @importFrom magrittr %>%
#' 
#' @examples
#' data <- POMA::st000336 # Example SummarizedExperiment object included in POMA
#' 
#' # T-test
#' results <- data %>%
#'   PomaImpute() %>% 
#'   PomaUnivariate() %>% 
#'   dplyr::select(feature, fold_change, pvalue)
#' 
#' results %>% 
#'   PomaVolcano(pval_cutoff = 0.05,
#'               log2fc_cutoff = NULL,
#'               labels = FALSE,
#'               x_label = "log2 (Fold Change)",
#'               y_label = "-log10 (P-value)")
#'               
#' # Limma
#' results <- data %>%
#'   PomaImpute() %>% 
#'   PomaNorm() %>% 
#'   PomaLimma(contrast = "DMD-Controls") %>% 
#'   dplyr::select(feature, log2FC, pvalue)
#' 
#' results %>% 
#'   PomaVolcano(pval_cutoff = 0.05,
#'               log2fc_cutoff = NULL,
#'               labels = FALSE,
#'               x_label = "log2 (Fold Change)",
#'               y_label = "-log10 (P-value)")
PomaVolcano <- function(data,
                        pval_cutoff = 0.05,
                        log2fc_cutoff = NULL,
                        labels = FALSE,
                        x_label = "log2 (Fold Change)",
                        y_label = "-log10 (P-value)") {



  to_volcano <- data %>% 
    as.data.frame() %>% 
    dplyr::select(feature = 1, logFC = 2, pvalue = 3)

  if (is.null(log2fc_cutoff)) {
    log2fc_cutoff <- quantile(abs(to_volcano$logFC), 0.75, na.rm = TRUE) 
  }
  
  plot_complete <- to_volcano %>% 
    ggplot2::ggplot(ggplot2::aes(logFC, -log10(pvalue), label = feature)) +
    ggplot2::geom_point(fill = "gray", size = 3, pch = 21, alpha = 0.8) +
    ggplot2::geom_point(data = to_volcano[to_volcano$pvalue < pval_cutoff & abs(to_volcano$logFC) > log2fc_cutoff ,], fill = "red", size = 3, pch = 21,) +
    {if(labels) ggrepel::geom_label_repel(data = to_volcano[to_volcano$pvalue < pval_cutoff & abs(to_volcano$logFC) > log2fc_cutoff ,], color = "black", size = 4)} +
    ggplot2::geom_vline(xintercept = c(-log2fc_cutoff, log2fc_cutoff), linetype = "dashed", color = "orange") +
    ggplot2::geom_hline(yintercept = -log10(pval_cutoff), linetype = "dashed", color = "orange") +
    ggplot2::labs(x = x_label,
                  y = y_label) +
    theme_poma()
  
  return(plot_complete)
}

