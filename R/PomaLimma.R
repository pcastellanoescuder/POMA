
#' Differential Expression Analysis Using `limma`
#'
#' @description `PomaLimma` uses the classical `limma` package to compute differential expression analysis.
#'
#' @param data A `SummarizedExperiment` object.
#' @param contrast Character. Indicates the comparison. For example, "Group1-Group2" or "control-intervention".
#' @param covs Character vector. Indicates the names of `colData` columns to be included as covariates. Default is NULL (no covariates). If not NULL, a limma model will be fitted using the specified covariates. Note: The order of the covariates is important and should be listed in increasing order of importance in the experimental design.
#' @param weights Logical. Indicates whether the limma model should estimate the relative quality weights for each group. See `?limma::arrayWeights()`.
#' @param replicates Character. Indicates the name of the `colData` column including the replicate factor label. Default is NULL (no replicates).
#' @param adjust Character. Indicates the multiple comparisons correction method. Options are: "fdr", "holm", "hochberg", "hommel", "bonferroni", "BH" and "BY".
#'
#' @export
#'
#' @return A `tibble` with the results.
#' @references Matthew E. Ritchie, Belinda Phipson, Di Wu, Yifang Hu, Charity W. Law, Wei Shi, Gordon K. Smyth, limma powers differential expression analyses for RNA-sequencing and microarray studies, Nucleic Acids Research, Volume 43, Issue 7, 20 April 2015, Page e47, https://doi.org/10.1093/nar/gkv007
#' @author Pol Castellano-Escuder
#'
#' @importFrom magrittr %>%
#' 
#' @examples 
#' data("st000284")
#' 
#' st000284 %>%
#'   PomaNorm() %>%
#'   PomaLimma(contrast = "Healthy-CRC", adjust = "fdr")
PomaLimma <- function(data,
                      contrast = NULL,
                      covs = NULL,
                      adjust = "fdr",
                      replicates = NULL,
                      weights = FALSE) {

  if (!is(data, "SummarizedExperiment")){
    stop("data is not a SummarizedExperiment object. \nSee POMA::PomaCreateObject or SummarizedExperiment::SummarizedExperiment")
  }
  if (ncol(SummarizedExperiment::colData(data)) == 0) {
    stop("metadata file required")
  }
  if (!is.factor(SummarizedExperiment::colData(data)[,1])) {
    stop("Grouping factor must be a factor (first column of the metadata file)")
  }
  if (is.null(contrast)) {
    stop("Specify a contrast")
  }
  if (!(adjust %in% c("fdr", "holm", "hochberg", "hommel", "bonferroni", "BH", "BY"))) {
    stop("Incorrect value for adjust argument")
  }

  main_factor <- as.factor(SummarizedExperiment::colData(data)[,1])
  to_limma <- SummarizedExperiment::assay(data)

  if (!is.null(covs)) {
    covariates <- SummarizedExperiment::colData(data) %>%
      as.data.frame() %>%
      dplyr::select(-1) %>%
      dplyr::select_at(dplyr::vars(dplyr::matches(covs)))
    
    form <- as.formula(noquote(paste("~ 0 + main_factor + ", paste0(colnames(covariates), collapse = " + ", sep = ""), sep = "")))
    
    design <- stats::model.matrix(form, data = covariates)
    colnames(design)[1:length(levels(main_factor))] <- levels(main_factor)
    
  } else {
    design <- stats::model.matrix( ~ 0 + main_factor)
    colnames(design) <- levels(main_factor)
  }
  
  if (!is.null(replicates)) {
    replicate <- SummarizedExperiment::colData(data) %>%
      as.data.frame() %>%
      dplyr::pull(replicates) %>% 
      factor()
    
    corfit <- limma::duplicateCorrelation(to_limma, design, block = replicate)
    lim_correlation <- corfit$consensus
  } else {
    replicate <- NULL
    lim_correlation <- NULL
  }
  
  cont_matrix <- limma::makeContrasts(contrasts = contrast, levels = design)
  
  if (weights) {
    array_weights <- limma::arrayWeights(to_limma, design = design)
    model <- limma::lmFit(to_limma, design, weights = array_weights, block = replicate, correlation = lim_correlation)
  } else {
    model <- limma::lmFit(to_limma, design, block = replicate, correlation = lim_correlation)
  }
  
  model <- limma::contrasts.fit(model, cont_matrix)
  
  modelstats <- limma::eBayes(model)
  result <- limma::topTable(modelstats, number = nrow(to_limma),
                            coef = contrast, sort.by = "p", adjust.method = adjust) %>% 
    tibble::rownames_to_column("feature") %>% 
    dplyr::rename(pvalue = P.Value, adj_pvalue = adj.P.Val) %>% 
    dplyr::as_tibble()
  
  return(result)
}

