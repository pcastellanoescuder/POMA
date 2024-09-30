
#' Univariate Statistical Test
#' 
#' @description `PomaUnivariate` performs parametric and non-parametric univariate statistical tests on a `SummarizedExperiment` object to compare groups or conditions. Available methods include T-test, ANOVA, ANCOVA, Mann Whitney U Test (Wilcoxon Rank Sum Test), and Kruskal-Wallis.
#'
#' @param data A `SummarizedExperiment` object.
#' @param method Character. The univariate statistical test to be performed. Available options include "ttest" (T-test), "anova" (analysis of variance), "mann" (Wilcoxon rank-sum test), and "kruskal" (Kruskal-Wallis test).
#' @param covs Character vector. Indicates the names of `colData` columns to be included as covariates. Default is NULL (no covariates). If not NULL, an ANCOVA model will be fitted using the specified covariates. Note: The order of the covariates is important and should be listed in increasing order of importance in the experimental design.
#' @param error Character vector. Indicates the name of a `colData` column to be included as an error term (e.g., replicates). Default is NULL (no error term).
#' @param paired Logical. Indicates if the data is paired or not. Default is FALSE.
#' @param var_equal Logical. Indicates if the data variances are assumed to be equal or not. Default is FALSE.
#' @param adjust Character. Multiple comparisons correction method to adjust p-values. Available options are: "fdr" (false discovery rate), "holm", "hochberg", "hommel", "bonferroni", "BH" (Benjamini-Hochberg), and "BY" (Benjamini-Yekutieli).
#' @param run_post_hoc Logical. Indicates if computing post-hoc tests or not. Setting this parameter to FALSE can save time for large datasets. 
#' 
#' @export
#'
#' @return A `tibble` for "ttest" and "mann". A `list` for "anova" and "kruskal".
#' @author Pol Castellano-Escuder
#'
#' @importFrom magrittr %>%
#' 
#' @examples 
#' # Two groups
#' ## Output columns: feature, fold_change, diff_means, pvalue, adj_pvalue, mean_xxx (group 1) mean_yyy (group 2), sd_xxx (group 1), sd_yyy (group 2)
#' data <- POMA::st000336 # Example SummarizedExperiment object included in POMA
#' 
#' ## Perform T-test
#' ttest_results <- st000336 %>% 
#'   PomaImpute() %>% 
#'   PomaUnivariate(method = "ttest",
#'                  paired = FALSE,
#'                  var_equal = FALSE,
#'                  adjust = "fdr")
#' 
#' ttest_results %>% 
#'   dplyr::slice(1:10)
#' 
#' ## Volcano plot
#' ttest_results %>% 
#'   dplyr::select(feature, fold_change, pvalue) %>% 
#'   PomaVolcano(labels = TRUE)
#' 
#' ## Boxplot of top features
#' data %>% 
#'   PomaBoxplots(x = "features", 
#'                outcome = "group", # factorial variable to group by (e.g., treatment, sex, etc)
#'                feature_name = ttest_results$feature[1:10])
#' 
#' ## Heatmap of top features
#' data[rownames(data) %in% ttest_results$feature[1:10]] %>% 
#'   PomaHeatmap(covs = c("group"), # covariates to plot (e.g., treatment, sex, etc)
#'               feature_names = TRUE)
#' 
#' ## Perform Mann-Whitney U test
#' mann_whitney_results <- st000336 %>% 
#'   PomaImpute() %>% 
#'   PomaUnivariate(method = "mann",
#'                  paired = FALSE,
#'                  var_equal = FALSE,
#'                  adjust = "fdr")
#' 
#' mann_whitney_results %>% 
#'   dplyr::slice(1:10)
#' 
#' ## Volcano plot
#' mann_whitney_results %>% 
#'   dplyr::select(feature, fold_change, pvalue) %>% 
#'   PomaVolcano(labels = TRUE)
#' 
#' ## Boxplot of top features
#' data %>% 
#'   PomaBoxplots(x = "features", 
#'                outcome = "group", # factorial variable to group by (e.g., treatment, sex, etc)
#'                feature_name = mann_whitney_results$feature[1:10])
#' 
#' ## Heatmap of top features
#' data[rownames(data) %in% mann_whitney_results$feature[1:10]] %>% 
#'   PomaHeatmap(covs = c("group"), # covariates to plot (e.g., treatment, sex, etc)
#'               feature_names = TRUE)
#' 
#' # More than 2 groups
#' ## Output is a list with objects `result` and `post_hoc_tests`
#' data <- POMA::st000284 # Example SummarizedExperiment object included in POMA
#' 
#' ## Perform Two-Way ANOVA
#' anova_results <- data %>% 
#'   PomaUnivariate(method = "anova",
#'                  covs = c("gender"),
#'                  error = NULL,
#'                  adjust = "fdr",
#'                  run_post_hoc = TRUE)
#' 
#' anova_results$result %>% 
#'   dplyr::slice(1:10)
#' 
#' anova_results$post_hoc_tests %>% 
#'   dplyr::slice(1:10)
#' 
#' ## Boxplot of top features
#' data %>% 
#'   PomaBoxplots(x = "features",
#'                outcome = "factors", # factorial variable to group by (e.g., treatment, sex, etc)
#'                feature_name = anova_results$result$feature[1:10])
#' 
#' ## Boxplot of top significant pairwise features (after posthoc test)
#' data %>% 
#'   PomaBoxplots(x = "features",
#'                outcome = "factors", # factorial variable to group by (e.g., treatment, sex, etc)
#'                feature_name = unique(anova_results$post_hoc_tests$feature)[1:10])
#' 
#' ## Heatmap of top features
#' data[rownames(data) %in% anova_results$result$feature[1:10]] %>% 
#'   PomaHeatmap(covs = c("factors"), # covariates to plot (e.g., treatment, sex, etc)
#'               feature_names = TRUE)
#' 
#' ## Perform Three-Way ANOVA
#' data %>% 
#'   PomaUnivariate(method = "anova", 
#'                  covs = c("gender", "smoking_condition"))
#' 
#' ## Perform ANCOVA with one numeric covariate and one factor covariate
#' data %>% 
#'   PomaUnivariate(method = "anova", 
#'                  covs = c("age_at_consent", "smoking_condition"))
#' 
#' # Perform Kruskal-Wallis test
#' data %>% 
#'   PomaUnivariate(method = "kruskal", 
#'                  adjust = "holm",
#'                  run_post_hoc = TRUE)
PomaUnivariate <- function(data,
                           method = "ttest",
                           covs = NULL,
                           error = NULL,
                           paired = FALSE,
                           var_equal = FALSE,
                           adjust = "fdr",
                           run_post_hoc = TRUE){

  if(!is(data, "SummarizedExperiment")){
    stop("data is not a SummarizedExperiment object. \nSee POMA::PomaCreateObject or SummarizedExperiment::SummarizedExperiment")
  }
  if (!(method %in% c("ttest", "anova", "mann", "kruskal"))) {
    stop("Incorrect value for method argument")
  }
  if (ncol(SummarizedExperiment::colData(data)) == 0) {
    stop("metadata file required")
  }
  if (!is.factor(SummarizedExperiment::colData(data)[,1])) {
    stop("Grouping factor must be a factor (first column of the metadata file)")
  }
  if (!(adjust %in% c("fdr", "holm", "hochberg", "hommel", "bonferroni", "BH", "BY"))) {
    stop("Incorrect value for adjust argument")
  }
  if (missing(method)) {
    message("method argument is empty. T-test will be used")
  }

  group_factor <- factor(SummarizedExperiment::colData(data)[,1])
  to_univariate <- t(SummarizedExperiment::assay(data))

  # group mean and SD
  group_means <- to_univariate %>%
    as.data.frame() %>% 
    dplyr::mutate(group = group_factor) %>%
    dplyr::group_by(group) %>%
    dplyr::summarise_all(list(~ mean(., na.rm = TRUE))) %>%
    tibble::remove_rownames() %>%
    tibble::column_to_rownames("group") %>%
    t() %>%
    as.data.frame() %>% 
    dplyr::rename_all(~ paste0("mean_", .))
  
  group_sd <- to_univariate %>%
    as.data.frame() %>% 
    dplyr::mutate(group = group_factor) %>%
    dplyr::group_by(group) %>%
    dplyr::summarise_all(list(~ sd(., na.rm = TRUE))) %>%
    tibble::remove_rownames() %>%
    tibble::column_to_rownames("group") %>%
    t() %>%
    as.data.frame() %>% 
    dplyr::rename_all(~ paste0("sd_", .))
  
  if (method == "ttest") {
    if (length(table(group_factor)[table(group_factor) != 0]) != 2) {
      stop("Grouping factor must have exactly 2 levels (first column of the metadata file)")
    }
    
    test_result <- data.frame(pvalue = apply(to_univariate, 2, function(x) {
      if (!paired) {
        res <- t.test(x ~ group_factor, data = data.frame(x, group_factor), 
                      na.rm = TRUE, 
                      alternative = "two.sided", 
                      var.equal = var_equal)
      } else {
        group1 <- x[group_factor == levels(group_factor)[1]]
        group2 <- x[group_factor == levels(group_factor)[2]]
        
        res <- t.test(group1, group2, 
                      na.rm = TRUE,
                      alternative = "two.sided", 
                      var.equal = var_equal, 
                      paired = TRUE)
      }
      res$p.value
    }))
    
    result <- test_result %>% 
      tibble::rownames_to_column("feature") %>%
      dplyr::mutate(adj_pvalue = p.adjust(pvalue, method = adjust)) %>%
      dplyr::bind_cols(group_means, group_sd) %>%
      dplyr::mutate(fold_change = as.numeric(group_means[,2] / group_means[,1]),
                    diff_means = as.numeric(group_means[,2] - group_means[,1])) %>%
      dplyr::select(feature, fold_change, diff_means, pvalue, adj_pvalue, dplyr::everything()) %>% 
      dplyr::arrange(pvalue) %>% 
      dplyr::as_tibble()

    return(result)
  }

  else if (method == "anova") {
    covariates <- SummarizedExperiment::colData(data) %>% 
      dplyr::as_tibble() %>% 
      dplyr::select(-1)
    
    if (is.null(covs)) {
      res_aov <- data.frame(pvalue = apply(to_univariate, 2, function(x) {anova(aov(x ~ group_factor))$"Pr(>F)"[1]})) %>%
        dplyr::mutate(adj_pvalue = p.adjust(pvalue, method = adjust)) %>% 
        dplyr::bind_cols(group_means, group_sd) %>% 
        tibble::rownames_to_column("feature") %>% 
        dplyr::select(feature, pvalue, adj_pvalue, dplyr::everything()) %>% 
        dplyr::as_tibble()

      # Post-hoc tests
      if (run_post_hoc) {
        post_hoc_tests <- list()
        for (i in seq_len(nrow(SummarizedExperiment::assay(data)))) {
          post_hoc_tests[[i]] <- dplyr::tibble(feature = rownames(SummarizedExperiment::assay(data))[i], 
                                               broom::tidy(TukeyHSD(aov(to_univariate[,i] ~ group_factor)))[,c(2, 7)])
        }
        
        post_hoc_tests <- dplyr::bind_rows(post_hoc_tests) %>% 
          dplyr::rename(adj_pvalue = adj.p.value) %>% 
          dplyr::arrange(adj_pvalue) 
      } else {
        post_hoc_tests <- NULL
      }
      
      return(list(result = res_aov, 
                  post_hoc_tests = post_hoc_tests))

    } else {
      if (is.null(error)) {
        covariates <- covariates %>%
          dplyr::select_at(dplyr::vars(dplyr::matches(covs)))
        
        model_names <- paste0(paste0(colnames(covariates), collapse = " * "), " * group_factor")
      } else {
        covariates <- covariates %>%
          dplyr::select_at(dplyr::vars(dplyr::matches(covs) | dplyr::matches(error)))
        
        model_names <- paste0(paste0(colnames(covariates)[colnames(covariates) != error], collapse = " * "), " * group_factor")
        model_names <- paste0(model_names, paste0(" + Error(", error,")"))
      }
      
      covariates_feat <- as.data.frame(cbind(to_univariate, covariates))
      result_cov <- vector(mode = "list", length = ncol(to_univariate))
      for(i in seq_len(ncol(to_univariate))) {
        result_cov[[i]] <- broom::tidy(aov(as.formula(paste(colnames(covariates_feat)[i], "~", model_names)), data = covariates_feat)) %>% 
          dplyr::filter(term != "Residuals" & !is.na(p.value)) %>% 
          dplyr::mutate(term = gsub("group_factor", names(SummarizedExperiment::colData(data))[1], term)) %>% 
          dplyr::select(term, pvalue = p.value) %>% 
          tidyr::pivot_wider(names_from = term, values_from = pvalue)
      }

      res_aov_cov <- result_cov %>% 
        dplyr::bind_rows() %>%
        dplyr::rename_all(~ paste0("pvalue_", .)) %>% 
        dplyr::mutate(dplyr::across(dplyr::starts_with("pvalue"), ~ p.adjust(., method = adjust), .names = "adj_pvalue_{.col}"),
                      feature = colnames(to_univariate)) %>% 
        dplyr::rename_at(dplyr::vars(dplyr::starts_with("adj_pvalue_pvalue_")), ~ gsub("adj_pvalue_pvalue", "adj_pvalue", .)) %>% 
        dplyr::bind_cols(group_means, group_sd) %>% 
        dplyr::select(feature, dplyr::everything()) %>% 
        dplyr::as_tibble()
      
      # Post-hoc tests
      if (run_post_hoc & is.null(error)) {
        post_hoc_tests <- list()
        for (i in seq_len(nrow(SummarizedExperiment::assay(data)))) {
          # Gereral Linear Hypotheses 
          posthoc <- multcomp::glht(
            aov(as.formula(paste(colnames(covariates_feat)[i], "~", model_names)),
                data = covariates_feat), 
            linfct = multcomp::mcp(group_factor = "Tukey"))
          summary_posthoc <- summary(posthoc)
          # Tukey
          # TukeyHSD(
          # aov(as.formula(paste(colnames(covariates_feat)[i], "~", model_names)),
          #     data = covariates_feat))$group_factor
          post_hoc_tests[[i]] <- dplyr::tibble(feature = rownames(SummarizedExperiment::assay(data))[i], 
                                               contrast = gsub(" ", "", names(summary_posthoc$test$coefficients)),
                                               pvalue = summary_posthoc$test$pvalues)
        }
        post_hoc_tests <- dplyr::bind_rows(post_hoc_tests) %>%
          dplyr::mutate(adj_pvalue = p.adjust(pvalue, method = "fdr")) %>% 
          dplyr::arrange(pvalue)
      } else {
        post_hoc_tests <- NULL
      }
      
      return(list(result = res_aov_cov, 
                  post_hoc_tests = post_hoc_tests))
    }
  }

  else if (method == "mann") {
    if (length(table(group_factor)[table(group_factor) != 0]) != 2) {
      stop("Grouping factor must have exactly 2 levels (first column of the metadata file)")
    }

    suppressWarnings({
      mann_result <- data.frame(pvalue = apply(to_univariate, 2, function(x) {
        if (!paired) {
          res <- wilcox.test(x ~ as.factor(group_factor), 
                             data = data.frame(x, as.factor(group_factor)))
        } else {
          group1 <- x[group_factor == levels(group_factor)[1]]
          group2 <- x[group_factor == levels(group_factor)[2]]
          
          res <- wilcox.test(group1, group2, 
                             paired = TRUE)
        }
        res$p.value
      }))
      
      result <- mann_result %>% 
        tibble::rownames_to_column("feature") %>%
        dplyr::mutate(adj_pvalue = p.adjust(pvalue, method = adjust)) %>%
        dplyr::bind_cols(group_means, group_sd) %>%
        dplyr::mutate(fold_change = as.numeric(group_means[,2]/group_means[,1]),
                      diff_means = as.numeric(group_means[,2] - group_means[,1])) %>% 
        dplyr::select(feature, fold_change, diff_means, pvalue, adj_pvalue, dplyr::everything()) %>% 
        dplyr::arrange(pvalue) %>% 
        dplyr::as_tibble()
    })
    
    return(result)
  }

  else if (method == "kruskal") {

    res_kruskal <- data.frame(pvalue = apply(to_univariate, 2, function(x){kruskal.test(x ~ as.factor(group_factor))$p.value})) %>%
      dplyr::mutate(adj_pvalue = p.adjust(pvalue, method = adjust),
                    kw_rank_sum = apply(to_univariate, 2, function(x){kruskal.test(x ~ as.factor(group_factor))$statistic})) %>%
      dplyr::bind_cols(group_means, group_sd) %>%
      tibble::rownames_to_column("feature") %>%
      dplyr::select(feature, kw_rank_sum, pvalue, adj_pvalue, dplyr::everything()) %>%
      dplyr::arrange(pvalue) %>% 
      dplyr::as_tibble()

    # Post-hoc tests
    if (run_post_hoc) {
      post_hoc_tests <- list()
      for (i in seq_len(nrow(SummarizedExperiment::assay(data)))) {
        post_hoc_tests[[i]] <- dplyr::tibble(feature = rownames(SummarizedExperiment::assay(data))[i],
                                             FSA::dunnTest(to_univariate[,i] ~ group_factor,
                                                           data = as.data.frame(to_univariate))$res
        )
      }
      
      post_hoc_tests <- dplyr::bind_rows(post_hoc_tests) %>%
        dplyr::select(feature, contrast = Comparison, adj_pvalue = P.adj) %>%
        dplyr::mutate(contrast = gsub(" ", "", contrast)) %>% 
        dplyr::arrange(adj_pvalue) 
    } else {
      post_hoc_tests <- NULL
    }
    
    return(list(result = res_kruskal, 
                post_hoc_tests = post_hoc_tests))
  }
}

