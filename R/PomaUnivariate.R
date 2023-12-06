
#' Univariate Statistical Test
#' 
#' @description `PomaUnivariate` performs parametric and non-parametric univariate statistical tests on a `SummarizedExperiment` object to compare groups or conditions. Available methods include T-test, ANOVA, ANCOVA, Mann Whitney U Test (Wilcoxon Rank Sum Test), and Kruskal-Wallis.
#'
#' @param data A `SummarizedExperiment` object.
#' @param method Character. The univariate statistical test to be performed. Available options include "ttest" (T-test), "anova" (analysis of variance), "mann" (Wilcoxon rank-sum test), and "kruskal" (Kruskal-Wallis test).
#' @param covs Character vector. Indicates the names of `colData` columns to be included as covariates. Default is NULL (no covariates). If not NULL, an ANCOVA model will be fitted using the specified covariates. Note: The order of the covariates is important and should be listed in increasing order of importance in the experimental design.
#' @param paired Logical. Indicates if the data is paired or not. Default is FALSE.
#' @param var_equal Logical. Indicates if the data variances are assumed to be equal or not. Default is FALSE.
#' @param adjust Character. Multiple comparisons correction method to adjust p-values. Available options are: "fdr" (false discovery rate), "holm", "hochberg", "hommel", "bonferroni", "BH" (Benjamini-Hochberg), and "BY" (Benjamini-Yekutieli).
#' @param run_post_hoc Logical. Indicates if computing post-hoc tests or not. Setting this parameter to FALSE can save time for large datasets. 
#' 
#' @export
#'
#' @return A `list` with the results.
#' @author Pol Castellano-Escuder
#'
#' @importFrom magrittr %>%
#' 
#' @examples 
#' data("st000336")
#' 
#' # Perform T-test
#' st000336 %>% 
#' PomaImpute() %>% 
#' PomaUnivariate(method = "ttest")
#' 
#' # Perform Mann-Whitney U test
#' st000336 %>% 
#' PomaImpute() %>% 
#' PomaUnivariate(method = "mann", paired = FALSE, adjust = "fdr")
#' 
#' data("st000284")
#' # Perform Two-Way ANOVA
#' st000284 %>% 
#' PomaUnivariate(method = "anova", covs = c("gender"))
#' 
#' # Perform Three-Way ANOVA
#' st000284 %>% 
#' PomaUnivariate(method = "anova", covs = c("gender", "smoking_condition"))
#' 
#' # Perform ANCOVA with one numeric covariate and one factor covariate
#' st000284 %>% 
#' PomaUnivariate(method = "anova", covs = c("age_at_consent", "smoking_condition"))
#' 
#' # Perform Kruskal-Wallis test
#' st000284 %>% 
#' PomaUnivariate(method = "kruskal", adjust = "holm")
PomaUnivariate <- function(data,
                           method = "ttest",
                           covs = NULL,
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

  group_factor <- SummarizedExperiment::colData(data)[,1]
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

    result <- data.frame(pvalue = apply(to_univariate, 2, function(x){t.test(x ~ group_factor, na.rm = TRUE, 
                                                                             alternative = "two.sided",
                                                                             var.equal = var_equal, 
                                                                             paired = paired)$p.value})) %>% 
      tibble::rownames_to_column("feature") %>%
      dplyr::mutate(adj_pvalue = p.adjust(pvalue, method = adjust)) %>%
      dplyr::bind_cols(group_means, group_sd) %>%
      dplyr::mutate(fold_change = as.numeric(round(group_means[,2] / group_means[,1], 3)),
                    diff_means = as.numeric(round(group_means[,2] - group_means[,1], 3))) %>%
      dplyr::select(feature, fold_change, diff_means, pvalue, adj_pvalue, dplyr::everything()) %>% 
      dplyr::arrange(pvalue) %>% 
      dplyr::as_tibble()

    return(list(result = result))
  }

  else if (method == "anova") {
    covariates <- SummarizedExperiment::colData(data) %>%
      as.data.frame() %>%
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
        for (i in 1:nrow(SummarizedExperiment::assay(data))) {
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
      covariates <- covariates %>%
        dplyr::select_at(dplyr::vars(dplyr::matches(covs)))

      model_names <- paste0(paste0(colnames(covariates), collapse = " * "), " * group_factor")
      covariates_feat <- as.data.frame(cbind(to_univariate, covariates))
      
      result_cov <- vector(mode = "list", length = ncol(to_univariate))
      for(i in 1:ncol(to_univariate)) {
        result_cov[[i]] <- broom::tidy(anova(aov(as.formula(paste(colnames(covariates_feat)[i], "~", model_names)),
                                                 data = covariates_feat))) %>% 
          dplyr::filter(term != "Residuals") %>% 
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
      if (run_post_hoc) {
        post_hoc_tests <- list()
        for (i in 1:nrow(SummarizedExperiment::assay(data))) {
          post_hoc_tests[[i]] <- dplyr::tibble(feature = rownames(SummarizedExperiment::assay(data))[i], 
                                               as.data.frame(TukeyHSD(
                                                 aov(as.formula(paste(colnames(covariates_feat)[1], "~", model_names)),
                                                     data = covariates_feat))$group_factor) %>% 
                                                 tibble::rownames_to_column("contrast") %>% 
                                                 dplyr::select(contrast, adj_pvalue = `p adj`)
          )
        }
        
        post_hoc_tests <- dplyr::bind_rows(post_hoc_tests) %>%
          dplyr::arrange(adj_pvalue)
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
      result <- data.frame(pvalue = apply(to_univariate, 2, function(x){wilcox.test(x ~ as.factor(group_factor),
                                                                                    paired = paired)$p.value})) %>% 
        tibble::rownames_to_column("feature") %>%
        dplyr::mutate(adj_pvalue = p.adjust(pvalue, method = adjust)) %>%
        dplyr::bind_cols(group_means, group_sd) %>%
        dplyr::mutate(fold_change = as.numeric(round(group_means[,2]/group_means[,1], 3)),
                      diff_means = as.numeric(round(group_means[,2] - group_means[,1], 3))) %>% 
        dplyr::select(feature, fold_change, diff_means, pvalue, adj_pvalue, dplyr::everything()) %>% 
        dplyr::arrange(pvalue) %>% 
        dplyr::as_tibble()
    })
    
    return(list(result = result))
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
      for (i in 1:nrow(SummarizedExperiment::assay(data))) {
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

