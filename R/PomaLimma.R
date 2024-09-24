
#' Differential Expression Analysis Using `limma`
#'
#' @description `PomaLimma` uses the classical `limma` package to compute differential expression analysis.
#'
#' @param data A `SummarizedExperiment` object.
#' @param contrast Character. Indicates the comparison. For example, "Group1-Group2" or "control-intervention".
#' @param outcome Character. Indicates the name of the `colData` column to be used as the outcome factor. Default is NULL (first factor variable in `colData`).
#' @param covs Character vector. Indicates the names of `colData` columns to be included as covariates. Default is NULL (no covariates). If not NULL, a limma model will be fitted using the specified covariates. Note: The order of the covariates is important and should be listed in increasing order of importance in the experimental design.
#' @param weights Logical. Indicates whether the limma model should estimate the relative quality weights for each group. See `?limma::arrayWeights()`.
#' @param block Character. Specifies the name of the `colData` factor column that includes the random effect variable to be considered (e.g., replicate). The default is NULL, indicating no random effect.
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
#' data <- POMA::st000284 %>% # Example SummarizedExperiment object included in POMA
#'   PomaNorm()
#' 
#' # Basic limma
#' limma_results <- data %>% 
#'   PomaLimma(contrast = "Healthy-CRC", 
#'             covs = NULL,
#'             adjust = "fdr",
#'             block = NULL)
#' 
#' limma_results %>% 
#'   dplyr::slice(1:10)
#' 
#' ## Volcano plot
#' limma_results %>% 
#'   dplyr::select(feature, log2FC, pvalue) %>% 
#'   PomaVolcano()
#' 
#' ## Boxplot of top features
#' data %>% 
#'   PomaBoxplots(x = "features", 
#'                outcome = "gender", # factorial variable to group by (e.g., treatment, sex, etc)
#'                feature_name = limma_results$feature[1:10])
#' 
#' ## Heatmap of top features
#' data[rownames(data) %in% limma_results$feature[1:10]] %>% 
#'   PomaHeatmap(covs = c("gender", "smoking_condition", "alcohol_consumption"), # covariates to plot (e.g., treatment, sex, etc)
#'               feature_names = TRUE)
#' 
#' # Basic limma on alternative outcome
#' SummarizedExperiment::colData(data)$gender <- factor(ifelse(SummarizedExperiment::colData(data)$gender == 0, "male", "female"))
#' data %>% 
#'   PomaLimma(contrast = "male-female", 
#'             outcome = "gender",
#'             covs = NULL,
#'             adjust = "fdr",
#'             block = NULL)
#' 
#' limma_results %>% 
#'   dplyr::slice(1:10)
#' 
#' ## Volcano plot
#' limma_results %>% 
#'   dplyr::select(feature, log2FC, pvalue) %>% 
#'   PomaVolcano()
#' 
#' ## Boxplot of top features
#' data %>% 
#'   PomaBoxplots(x = "features", 
#'                outcome = "gender", # factorial variable to group by (e.g., treatment, sex, etc)
#'                feature_name = limma_results$feature[1:10])
#' 
#' ## Heatmap of top features
#' data[rownames(data) %in% limma_results$feature[1:10]] %>% 
#'   PomaHeatmap(covs = c("gender", "smoking_condition", "alcohol_consumption"), # covariates to plot (e.g., treatment, sex, etc)
#'               feature_names = TRUE)
#' 
#' # limma with one covariate
#' data %>% 
#'   PomaLimma(contrast = "Healthy-CRC", 
#'             covs = "gender",
#'             adjust = "fdr",
#'             block = NULL)
#' 
#' limma_results %>% 
#'   dplyr::slice(1:10)
#' 
#' ## Volcano plot
#' limma_results %>% 
#'   dplyr::select(feature, log2FC, pvalue) %>% 
#'   PomaVolcano()
#' 
#' ## Boxplot of top features
#' data %>% 
#'   PomaBoxplots(x = "features", 
#'                outcome = "gender", # factorial variable to group by (e.g., treatment, sex, etc)
#'                feature_name = limma_results$feature[1:10])
#' 
#' ## Heatmap of top features
#' data[rownames(data) %in% limma_results$feature[1:10]] %>% 
#'   PomaHeatmap(covs = c("gender", "smoking_condition", "alcohol_consumption"), # covariates to plot (e.g., treatment, sex, etc)
#'               feature_names = TRUE)
#' 
#' # limma with two covariates
#' data %>% 
#'   PomaLimma(contrast = "Healthy-CRC", 
#'             covs = c("gender", "age_at_consent"),
#'             adjust = "fdr",
#'             block = NULL)
#' 
#' limma_results %>% 
#'   dplyr::slice(1:10)
#' 
#' ## Volcano plot
#' limma_results %>% 
#'   dplyr::select(feature, log2FC, pvalue) %>% 
#'   PomaVolcano()
#' 
#' ## Boxplot of top features
#' data %>% 
#'   PomaBoxplots(x = "features", 
#'                outcome = "gender", # factorial variable to group by (e.g., treatment, sex, etc)
#'                feature_name = limma_results$feature[1:10])
#' 
#' ## Heatmap of top features
#' data[rownames(data) %in% limma_results$feature[1:10]] %>% 
#'   PomaHeatmap(covs = c("gender", "smoking_condition", "alcohol_consumption"), # covariates to plot (e.g., treatment, sex, etc)
#'               feature_names = TRUE)
#' 
#' # limma with replicates
#' # data %>% 
#' #   PomaLimma(contrast = "Healthy-CRC", 
#' #             covs = NULL,
#' #             adjust = "fdr",
#' #             block = "replicate")
PomaLimma <- function(data,
                      contrast = NULL,
                      outcome = NULL,
                      covs = NULL,
                      adjust = "fdr",
                      block = NULL,
                      weights = FALSE) {

  if (!is(data, "SummarizedExperiment")){
    stop("data is not a SummarizedExperiment object. \nSee POMA::PomaCreateObject or SummarizedExperiment::SummarizedExperiment")
  }
  if (ncol(SummarizedExperiment::colData(data)) == 0) {
    stop("metadata file required")
  }
  if (is.null(contrast)) {
    stop("Specify a contrast")
  }
  if (!(adjust %in% c("fdr", "holm", "hochberg", "hommel", "bonferroni", "BH", "BY"))) {
    stop("Incorrect value for adjust argument")
  }

  if (is.null(outcome)) {
    grouping_factor <- SummarizedExperiment::colData(data)[,1]
  } else {
    grouping_factor <- SummarizedExperiment::colData(data) %>%
      as.data.frame() %>%
      dplyr::pull(outcome)
  }
  
  if (!is.factor(grouping_factor)) {
    stop("Outcome (dependent variable) must be a factor")
  }
  
  grouping_factor <- factor(grouping_factor)
  to_limma <- SummarizedExperiment::assay(data)

  if (!is.null(covs)) {
    covariates <- SummarizedExperiment::colData(data) %>%
      as.data.frame() %>%
      dplyr::select_at(dplyr::vars(dplyr::matches(covs)))
    
    form <- as.formula(noquote(paste("~ 0 + grouping_factor + ", paste0(colnames(covariates), collapse = " + ", sep = ""), sep = "")))
    
    design <- stats::model.matrix(form, data = covariates)
    colnames(design)[1:length(levels(grouping_factor))] <- levels(grouping_factor)
    
  } else {
    design <- stats::model.matrix( ~ 0 + grouping_factor)
    colnames(design) <- levels(grouping_factor)
  }
  
  if (!is.null(block)) {
    replicate <- SummarizedExperiment::colData(data) %>%
      as.data.frame() %>%
      dplyr::pull(block) %>% 
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
    dplyr::rename(log2FC = logFC, pvalue = P.Value, adj_pvalue = adj.P.Val) %>% 
    dplyr::as_tibble()
  
  return(result)
}

