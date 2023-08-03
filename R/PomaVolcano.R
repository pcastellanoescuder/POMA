
#' Volcano Plot
#'
#' @description `PomaVolcano` generates a volcano plot. Data should not contain negative values.  
#'
#' @param data A `SummarizedExperiment` object.
#' @param method Character. Indicates the statistical method to use. Options are "ttest", "mann", "limma", and "DESeq".
#' @param pval Character. Indicates the pvalue type to generate the volcano plot. Options are: "raw" and "adjusted".
#' @param pval_cutoff Numeric. Indicated the pvalue cutoff (horizontal line).
#' @param adjust Character, Indicates the multiple correction method. Options are: "fdr", "holm", "hochberg", "hommel", "bonferroni", "BH" and "BY".
#' @param log2fc_cutoff Numeric. Indicates the log2 fold change cutoff (vertical lines).
#' @param labels Logical. Indicates if significant labels should be plotted. Defaul is FALSE.
#' @param paired Logical. Indicates if the data is paired or not. Default is FALSE.
#' @param var_equal Logical. Indicates if the data variances are assumed to be equal or not. Default is FALSE.
#' 
#' @export
#'
#' @return A `ggplot` object.
#' @author Pol Castellano-Escuder
#'
#' @importFrom magrittr %>%
#' 
#' @examples 
#' data("st000336")
#' 
#' st000336 %>% 
#'   PomaImpute() %>%
#'   PomaVolcano()
#'   
#' st000336 %>% 
#'   PomaImpute() %>%
#'   PomaVolcano(labels = TRUE)
#'   
#' st000336 %>% 
#'   PomaImpute() %>%
#'   PomaVolcano(labels = TRUE, pval = "adjusted")
PomaVolcano <- function(data,
                        method = "ttest",
                        pval = "raw",
                        pval_cutoff = 0.05,
                        adjust = "fdr",
                        log2fc_cutoff = NULL,
                        labels = FALSE,
                        paired = FALSE,
                        var_equal = FALSE,
                        ...) {

  if(!is(data, "SummarizedExperiment")){
    stop("data is not a SummarizedExperiment object. \nSee POMA::PomaCreateObject or SummarizedExperiment::SummarizedExperiment")
  }
  if (length(table(SummarizedExperiment::colData(data)[,1])) != 2) {
    stop("Grouping factor must have exactly 2 levels (first column of the metadata file)")
  }
  if (!(pval %in% c("raw", "adjusted"))) {
    stop('Incorrect value for pval argument. Options are: "raw" and "adjusted"')
  }
  if (!(adjust %in% c("fdr", "holm", "hochberg", "hommel", "bonferroni", "BH", "BY"))) {
    stop("Incorrect value for adjust argument.")
  }

  if (method == "ttest") {
    volcano_res <- data %>% 
      PomaUnivariate(method = "ttest", adjust = adjust, paired = paired, var_equal = var_equal) %>%
      dplyr::mutate(logFC = log2(fold_change))
  }
  else if (method == "mann") {
    volcano_res <- data %>% 
      PomaUnivariate(method = "mann", adjust = adjust, paired = paired, var_equal = var_equal) %>%
      dplyr::mutate(logFC = log2(fold_change))
  }
  else if (method == "limma") {
    contrast <- paste0(names(table(SummarizedExperiment::colData(data)[,1])), collapse = "-")
    
    volcano_res <- data %>% 
      PomaLimma(adjust = adjust, contrast = contrast)
  }
  else if (method == "DESeq") {
    volcano_res <- data %>%
      PomaDESeq(adjust = adjust)
  }
  
  if (pval == "raw") {
    volcano_res <- volcano_res %>% 
      dplyr::select(feature, logFC, pvalue = pvalue)
  } else {
    volcano_res <- volcano_res %>% 
      dplyr::select(feature, logFC, pvalue = adj_pvalue)
  }

  if (is.null(log2fc_cutoff)) {
    log2fc_cutoff <- quantile(abs(volcano_res$logFC), 0.75) 
  }
  
  plot_complete <- volcano_res %>% 
    ggplot2::ggplot(ggplot2::aes(logFC, -log10(pvalue), label = feature)) +
    ggplot2::geom_point(fill = "gray", size = 3, pch = 21, alpha = 0.6) +
    ggplot2::geom_point(data = volcano_res[volcano_res$pvalue < pval_cutoff & abs(volcano_res$logFC) > log2fc_cutoff ,], fill = "red", size = 3, pch = 21,) +
    {if(labels) ggrepel::geom_label_repel(data = volcano_res[volcano_res$pvalue < pval_cutoff & abs(volcano_res$logFC) > log2fc_cutoff ,], color = "black", size = 4)} +
    ggplot2::geom_vline(xintercept = c(-log2fc_cutoff, log2fc_cutoff), linetype = "dashed", color = "orange") +
    ggplot2::geom_hline(yintercept = -log10(pval_cutoff), linetype = "dashed", color = "orange") +
    ggplot2::labs(x = "log2 (Fold Change)",
                  y = ifelse(pval == "raw", "-log10 (P-value)", "-log10 (Adjusted P-value)")) +
    theme_poma()
  
  return(plot_complete)
  
}

