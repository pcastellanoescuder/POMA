
#' Differential Expression Analysis Based on the Negative Binomial Distribution
#'
#' @description `PomaDESeq` estimates variance-mean dependence in count data from high-throughput sequencing assays and test for differential expression based on a model using the negative binomial distribution.
#'
#' @param data A `SummarizedExperiment` object.
#' @param contrast Character. Indicates the comparison. For example, "Group1-Group2" or "control-intervention".
#' @param outcome Character. Indicates the name of the `colData` column to be used as the outcome factor. Default is NULL (first factor variable in `colData`).
#' @param covs Character vector. Indicates the names of `colData` columns to be included as covariates. Default is NULL (no covariates). If not NULL, a limma model will be fitted using the specified covariates. Note: The order of the covariates is important and should be listed in increasing order of importance in the experimental design.
#' @param adjust Character. Indicates the multiple comparisons correction method. Options are: "fdr", "holm", "hochberg", "hommel", "bonferroni", "BH" and "BY".
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
#' #library("airway")
#' #data("airway")
#' #se <- airway
#' #
#' ## Classic DESeq2
#' #se %>% 
#' #  PomaDESeq(contrast = NULL,
#' #            outcome = "dex",
#' #            covs = NULL,
#' #            adjust = "fdr")
#' #
#' ## DESeq2 with covariates
#' #se %>% 
#' #  PomaDESeq(contrast = NULL,
#' #            outcome = "dex",
#' #            covs = "cell",
#' #            adjust = "fdr")
#' #
#' ## DESeq2 with covariates and batch
#' #se %>% 
#' #  PomaDESeq(contrast = NULL,
#' #            outcome = "dex",
#' #            covs = c("batch", "cell"),
#' #            adjust = "fdr")
PomaDESeq <- function(data,
                      contrast = NULL,
                      outcome = NULL,
                      covs = NULL,
                      adjust = "fdr") {
  
  if (!is(data, "SummarizedExperiment")){
    stop("data is not a SummarizedExperiment object. \nSee POMA::PomaCreateObject or SummarizedExperiment::SummarizedExperiment")
  }
  if (ncol(SummarizedExperiment::colData(data)) == 0) {
    stop("metadata file required")
  }
  if (!(adjust %in% c("fdr", "holm", "hochberg", "hommel", "bonferroni", "BH", "BY"))) {
    stop("Incorrect value for adjust argument")
  }
  
  if (is.null(outcome)) {
    outcome <- colnames(SummarizedExperiment::colData(data))[1]
  }
  
  counts <- SummarizedExperiment::assay(data)
  coldata <- SummarizedExperiment::colData(data) %>% 
    as.data.frame() %>% 
    dplyr::select(dplyr::all_of(c(outcome, covs)))
  
  # Combine covariates and outcome into a single formula
  formula_expr <- purrr::reduce(rlang::syms(c(covs, outcome)), function(x, y) rlang::expr(!!x + !!y))
  
  deseq_object <- DESeq2::DESeqDataSetFromMatrix(countData = counts,
                                                 colData = coldata,
                                                 design = formula(rlang::expr(~ !!formula_expr))
                                                 )
  
  if (is.null(contrast)) {
    contrast <- coldata %>% 
      dplyr::pull(outcome) %>% 
      table() %>% 
      names() %>% 
      paste0(collapse = "-")
  }
  
  result <- deseq_object %>% 
    DESeq2::DESeq() %>% 
    DESeq2::results(pAdjustMethod = adjust, 
                    contrast = c(outcome, 
                                 unlist(strsplit(contrast, split = "-"))[1], 
                                 unlist(strsplit(contrast, split = "-"))[2])) %>% 
    dplyr::as_tibble(rownames = "feature") %>% 
    dplyr::rename(adj_pvalue = padj) %>% 
    dplyr::arrange(adj_pvalue)
  
  return(result)
}

