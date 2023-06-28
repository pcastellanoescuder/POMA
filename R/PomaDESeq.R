
#' Differential Expression Analysis Based on the Negative Binomial Distribution
#'
#' @description `PomaDESeq` is a `DESeq2` package wrapper to estimate variance-mean dependence in count data from high-throughput sequencing assays and test for differential expression based on a model using the negative binomial distribution.
#'
#' @param data A `SummarizedExperiment` object.
#' @param adjust Multiple comparisons correction method. Options are: "fdr", "holm", "hochberg", "hommel", "bonferroni", "BH", and "BY".
#' 
#' @export
#'
#' @return A `tibble` with the results.
#' @references Love, M.I., Huber, W., Anders, S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2 Genome Biology 15(12):550 (2014)
#' @author Pol Castellano-Escuder
#'
#' @importFrom magrittr %>%
#' 
#' @examples
#' data("st000284")
#' 
#' st000284_sub <- st000284[, st000284@colData$factors %in% c("CRC", "Healthy")] # select two groups
#' 
#' SummarizedExperiment::assay(st000284_sub) <- floor(SummarizedExperiment::assay(st000284_sub)) # convert all values to integers
#' 
#' st000284_sub %>% 
#'   PomaDESeq()
PomaDESeq <- function(data,
                      adjust = "fdr") {
  
  if(!is(data, "SummarizedExperiment")){
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
  
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = counts,
                                        colData = coldata,
                                        design = ~ condition
                                        )

  res_df <- dds %>% 
    DESeq2::DESeq() %>% 
    DESeq2::results(pAdjustMethod = adjust) %>% 
    dplyr::as_tibble(rownames = "feature") %>% 
    dplyr::rename(adj_pvalue = padj) %>% 
    dplyr::arrange(adj_pvalue)
  
  return(res_df)
}

