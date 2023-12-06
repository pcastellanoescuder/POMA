
#' Differential Expression Analysis Based on the Negative Binomial Distribution
#'
#' @description `PomaDESeq` estimates variance-mean dependence in count data from high-throughput sequencing assays and test for differential expression based on a model using the negative binomial distribution.
#'
#' @param data A `SummarizedExperiment` object.
#' @param adjust Character. Multiple comparisons correction method to adjust p-values. Available options are: "fdr" (false discovery rate), "holm", "hochberg", "hommel", "bonferroni", "BH" (Benjamini-Hochberg), and "BY" (Benjamini-Yekutieli).
#' 
#' @export
#'
#' @return A `tibble` with the results.
#' @references Love, M.I., Huber, W., Anders, S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2 Genome Biology 15(12):550 (2014)
#' @author Pol Castellano-Escuder
#'
#' @importFrom magrittr %>%
PomaDESeq <- function(data,
                      adjust = "fdr") {
  
  if (!is(data, "SummarizedExperiment")){
    stop("data is not a SummarizedExperiment object. \nSee POMA::PomaCreateObject or SummarizedExperiment::SummarizedExperiment")
  }
  if (ncol(SummarizedExperiment::colData(data)) == 0) {
    stop("metadata file required")
  }
  if (!is.factor(SummarizedExperiment::colData(data)[,1])) {
    stop("Grouping factor must be a factor (first column of the metadata file)")
  }
  
  counts <- SummarizedExperiment::assay(data)
  coldata <- SummarizedExperiment::colData(data) %>% 
    as.data.frame() %>% 
    tibble::rownames_to_column("sample_id") %>% 
    dplyr::rename(condition = 2)
  
  deseq_object <- DESeq2::DESeqDataSetFromMatrix(countData = counts,
                                                 colData = coldata,
                                                 design = ~ condition
                                                 )

  res_df <- deseq_object %>% 
    DESeq2::DESeq() %>% 
    DESeq2::results(pAdjustMethod = adjust) %>% 
    dplyr::as_tibble(rownames = "feature") %>% 
    dplyr::rename(adj_pvalue = padj) %>% 
    dplyr::arrange(adj_pvalue)
  
  return(res_df)
}

