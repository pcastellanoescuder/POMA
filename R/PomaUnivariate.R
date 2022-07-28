
#' Univariate Statistical Methods for Mass Spectrometry Data
#'
#' @description PomaUnivariate() allows users to perform different univariate statistical analysis on MS data.
#'
#' @param data A SummarizedExperiment object.
#' @param covariates Logical. If it's set to `TRUE` all metadata variables stored in `colData` will be used as covariables. Default = FALSE.
#' @param covs Character vector indicating the name of `colData` columns that will be included as covariates. Default is NULL (all variables).
#' @param method Univariate statistical method. Options are: "ttest", "anova", "mann" and "kruskal".
#' @param paired Logical that indicates if the data is paired or not.
#' @param var_equal Logical that indicates if the data variance is equal or not.
#' @param adjust Multiple comparisons correction method. Options are: "fdr", "holm", "hochberg", "hommel", "bonferroni", "BH" and "BY".
#'
#' @export
#'
#' @return A tibble with results.
#' @author Pol Castellano-Escuder
#'
#' @importFrom tibble column_to_rownames rownames_to_column remove_rownames
#' @importFrom dplyr select mutate mutate_all rename_at filter bind_cols bind_rows summarise_all group_by as_tibble everything rename_all across starts_with vars
#' @importFrom magrittr %>%
#' @importFrom SummarizedExperiment assay colData
#' 
#' @examples 
#' data("st000336")
#' data("st000284")
#' 
#' # ttest
#' st000336 %>% 
#'   PomaImpute() %>%
#'   PomaNorm() %>%
#'   PomaOutliers() %>%
#'   PomaUnivariate(method = "ttest")
#' 
#' # ANOVA
#' st000284 %>% 
#'   PomaImpute() %>%
#'   PomaNorm() %>%
#'   PomaOutliers() %>%
#'   PomaUnivariate(method = "anova")
PomaUnivariate <- function(data,
                           covariates = FALSE,
                           covs = NULL,
                           method = "ttest",
                           paired = FALSE,
                           var_equal = FALSE,
                           adjust = "fdr"){

  if (missing(data)) {
    stop("data argument is empty!")
  }
  if(!is(data[1], "SummarizedExperiment")){
    stop("data is not a SummarizedExperiment object. \nSee POMA::PomaSummarizedExperiment or SummarizedExperiment::SummarizedExperiment")
  }
  if (missing(method)) {
    stop("Select a method!")
  }
  if (!(method %in% c("ttest", "anova", "mann", "kruskal"))) {
    stop("Incorrect value for method argument!")
  }
  if (!(adjust %in% c("fdr", "holm", "hochberg", "hommel", "bonferroni", "BH", "BY"))) {
    stop("Incorrect value for adjust argument!")
  }

  if(covariates & ncol(SummarizedExperiment::colData(data)) == 1){
    stop("Seems there aren't covariates in your data...")
  }

  Group <- as.factor(SummarizedExperiment::colData(data)[,1])
  e <- t(SummarizedExperiment::assay(data))

  ## group means and sd
  
  group_means <- e %>%
    as.data.frame() %>% 
    mutate(group = Group) %>%
    group_by(group) %>%
    summarise_all(list(~ mean(., na.rm = TRUE))) %>%
    remove_rownames() %>%
    column_to_rownames("group") %>%
    t() %>%
    as.data.frame() %>% 
    rename_all(~ paste0("mean_", .))
  
  group_sd <- e %>%
    as.data.frame() %>% 
    mutate(group = Group) %>%
    group_by(group) %>%
    summarise_all(list(~ sd(., na.rm = TRUE))) %>%
    remove_rownames() %>%
    column_to_rownames("group") %>%
    t() %>%
    as.data.frame() %>% 
    rename_all(~ paste0("sd_", .))
  
  ##
  
  if(method == "ttest"){

    stat_ttest <- function(x){t.test(x ~ Group, na.rm = TRUE, alternative = "two.sided",
                                     var.equal = var_equal, paired = paired)$p.value}

    res_ttest <- data.frame(pvalue = apply(FUN = stat_ttest, MARGIN = 2, X = e)) %>% 
      rownames_to_column("feature") %>%
      mutate(pvalueAdj = p.adjust(pvalue, method = adjust)) %>%
      bind_cols(group_means, group_sd) %>%
      mutate(FC = as.numeric(round(group_means[,2]/group_means[,1], 3)),
             diff_means = as.numeric(round(group_means[,2] - group_means[,1], 3))) %>%
      dplyr::select(feature, FC, diff_means, pvalue, pvalueAdj, everything()) %>% 
      dplyr::as_tibble()

    return(res_ttest)
  }

  else if(method == "anova"){

    if(!covariates){

      stat_aov <- function(x){anova(aov(x ~ Group))$"Pr(>F)"[1]}
      
      res_aov <- data.frame(pvalue = apply(FUN = stat_aov, MARGIN = 2, X = e)) %>%
        mutate(pvalueAdj = p.adjust(pvalue, method = adjust)) %>% 
        bind_cols(group_means, group_sd) %>% 
        rownames_to_column("feature") %>% 
        dplyr::select(feature, pvalue, pvalueAdj, everything()) %>% 
        dplyr::as_tibble()

      return(res_aov)

    }
    else {
      
      if(is.null(covs)){
        covariates <- colData(data) %>%
          as.data.frame() %>%
          dplyr::select(-1) %>% 
          dplyr::mutate_all(as.numeric)
      } 
      else {
        covariates <- colData(data) %>%
          as.data.frame() %>%
          dplyr::select(-1) %>% 
          dplyr::select_at(vars(matches(covs))) %>% 
          dplyr::mutate_all(as.numeric)
      }

      model_names <- paste0(paste0(colnames(covariates), collapse = " + "), " + Group")
      covariates_feat <- as.data.frame(cbind(e, covariates))
      
      result_cov <- vector(mode = "list", length = ncol(e))
      for(i in 1:ncol(e)) {
        result_cov[[i]] <- anova(aov(as.formula(paste(colnames(covariates_feat)[i], "~", 
                                                      model_names)),
                                     data = covariates_feat))$"Pr(>F)"[1:(ncol(covariates)+1)]
        
        names(result_cov[[i]]) <- c(colnames(covariates), colnames(SummarizedExperiment::colData(data))[1])
      }

      res_aov_cov <- bind_rows(result_cov) %>%
        dplyr::rename_all(~ paste0("pvalue_", .)) %>% 
        dplyr::mutate(across(starts_with("pvalue"), ~ p.adjust(., method = adjust), .names = "pvalueAdj_{.col}"),
                      feature = colnames(e)) %>% 
        dplyr::rename_at(vars(starts_with("pvalueAdj_pvalue_")), ~ gsub("pvalueAdj_pvalue", "pvalueAdj", .)) %>% 
        bind_cols(group_means, group_sd) %>% 
        dplyr::select(feature, everything()) %>% 
        dplyr::as_tibble()

      return(res_aov_cov)

    }
  }

  else if(method == "mann"){

    res_mann <- data.frame(pvalue = apply(e, 2, function(x){wilcox.test(x ~ as.factor(Group),
                                                                        paired = paired)$p.value})) %>% 
      rownames_to_column("feature") %>%
      mutate(pvalueAdj = p.adjust(pvalue, method = adjust)) %>%
      bind_cols(group_means, group_sd) %>%
      mutate(FC = as.numeric(round(group_means[,2]/group_means[,1], 3)),
             diff_means = as.numeric(round(group_means[,2] - group_means[,1], 3))) %>% 
      dplyr::select(feature, FC, diff_means, pvalue, pvalueAdj, everything()) %>% 
      dplyr::as_tibble()
    
    return(res_mann)
  }

  else if (method == "kruskal"){

    res_kruskal <- data.frame(pvalue = apply(e, 2, function(x){kruskal.test(x ~ as.factor(Group))$p.value})) %>%
      mutate(pvalueAdj = p.adjust(pvalue, method = adjust),
             kw_rank_sum = apply(e, 2, function(x){kruskal.test(x ~ as.factor(Group))$statistic})) %>% 
      bind_cols(group_means, group_sd) %>% 
      rownames_to_column("feature") %>%
      dplyr::select(feature, kw_rank_sum, pvalue, pvalueAdj, everything()) %>% 
      dplyr::as_tibble()
    
    return(res_kruskal)
  }

}

