
#' Univariate Statistical Methods for Mass Spectrometry Data
#'
#' @description PomaUnivariate() allows users to perform different univariate statistical analysis on MS data.
#'
#' @param data A MSnSet object. First `pData` column must be the subject group/type.
#' @param covariates Logical. If it's set to `TRUE` all metadata variables stored in `pData` will be used as covariables. Default = FALSE.
#' @param method Univariate statistical method. Options are: "ttest", "anova", "mann" and "kruskal".
#' @param paired Logical that indicates if the data is paired or not.
#' @param var_equal Logical that indicates if the data variance is equal or not.
#' @param adjust Multiple comparisons correction method. Options are: "fdr", "holm", "hochberg", "hommel", "bonferroni", "BH" and "BY".
#'
#' @export
#'
#' @return A data frame with results.
#' @author Pol Castellano-Escuder
#'
#' @importFrom tibble column_to_rownames rownames_to_column
#' @importFrom dplyr select mutate filter bind_cols bind_rows 
#' @importFrom magrittr %>%
#' @importFrom crayon red
#' @importFrom clisymbols symbol
#' @importFrom Biobase varLabels pData exprs
PomaUnivariate <- function(data,
                           covariates = FALSE,
                           method = "ttest",
                           paired = FALSE,
                           var_equal = FALSE,
                           adjust = "fdr"){

  if (missing(data)) {
    stop(crayon::red(clisymbols::symbol$cross, "data argument is empty!"))
  }
  if(!(class(data) == "MSnSet")){
    stop(paste0(crayon::red(clisymbols::symbol$cross, "data is not a MSnSet object."), 
                " \nSee POMA::PomaMSnSetClass or MSnbase::MSnSet"))
  }
  if (missing(method)) {
    stop(crayon::red(clisymbols::symbol$cross, "Select a method!"))
  }
  if (!(method %in% c("ttest", "anova", "mann", "kruskal"))) {
    stop(crayon::red(clisymbols::symbol$cross, "Incorrect value for method argument!"))
  }
  if (!(adjust %in% c("fdr", "holm", "hochberg", "hommel", "bonferroni", "BH", "BY"))) {
    stop(crayon::red(clisymbols::symbol$cross, "Incorrect value for adjust argument!"))
  }

  if(isTRUE(covariates) & ncol(pData(data)) == 1){
    stop(crayon::red(clisymbols::symbol$cross, "Seems that your data don't have covariates..."))
  }

  Biobase::varLabels(data)[1] <- "Group"
  Group <- as.factor(Biobase::pData(data)$Group)
  e <- t(Biobase::exprs(data))

  ## calcule means
  
  group_means <- e %>% 
    as.data.frame() %>% 
    mutate(group = Group)
  
  suppressWarnings({
    group_means <- data.frame(aggregate(group_means, by = list(group_means$group), FUN = mean)) %>% 
      column_to_rownames("Group.1") %>% 
      t() %>% 
      as.data.frame() %>% 
      round(2) %>% 
      filter(rownames(.) != "group")
    
    colnames(group_means) <- paste0("mean_", colnames(group_means))
  })
  
  ##
  
  if(method == "ttest"){

    stat <- function(x){t.test(x ~ Group, na.rm = TRUE, alternative = c("two.sided"),
                               var.equal = var_equal, paired = paired)$p.value}

    p <- data.frame(pvalue = apply(FUN = stat, MARGIN = 2, X = e))

    p <- p %>%
      rownames_to_column("ID") %>%
      mutate(pvalueAdj = p.adjust(pvalue, method = adjust)) %>%
      column_to_rownames("ID")
    
    p <- bind_cols(group_means, p) %>%
      rownames_to_column("ID") %>%
      mutate(Fold_Change_Ratio = as.numeric(round(group_means[,2]/group_means[,1], 3)),
             Difference_Of_Means = as.numeric(round(group_means[,2] - group_means[,1], 3))) %>%
      column_to_rownames("ID") %>%
      select(1,2,5,6,3,4)

    return(p)
  }

  else if(method == "anova"){

    if(!isTRUE(covariates)){

      stat2 <- function(x){anova(aov(x ~ Group))$"Pr(>F)"[1]}
      p2 <- data.frame(pvalue = apply(FUN = stat2, MARGIN = 2, X = e))

      p2 <- p2 %>%
        mutate(pvalueAdj = p.adjust(pvalue, method = adjust))
      
      p2 <- bind_cols(group_means, p2)

      return(p2)

    }
    else{

      covariate_uni <- pData(data)[, 2:ncol(pData(data))]
      covariate_uni <- sapply(covariate_uni, as.numeric)

      model_names <- paste0(paste0(colnames(covariate_uni), collapse = " + "), " + Group")

      LenCov <- ncol(covariate_uni)

      covariate_uni <- as.data.frame(cbind(e, covariate_uni))

      n <- ncol(covariate_uni) - LenCov
      result <- vector(mode = "list", length = n)
      for(i in 1:n) {
        result[[i]] <- data.frame(pvalue = anova(aov(as.formula(paste(colnames(covariate_uni)[i], "~", model_names)),
                                            data = covariate_uni))$"Pr(>F)"[LenCov+1])
      }

      p3 <- bind_rows(result)
      rownames(p3) <- colnames(e)

      p3 <- p3 %>%
        mutate(pvalueAdj = p.adjust(pvalue, method = adjust))

      p3 <- bind_cols(group_means, p3)

      return(p3)

    }
  }

  else if(method == "mann"){

    non_param_mann <- data.frame(pvalue = apply(e, 2,
                                                 function(x){wilcox.test(x ~ as.factor(Group),
                                                                         paired = paired)$p.value}))

    non_param_mann <- non_param_mann %>%
      rownames_to_column("ID") %>%
      mutate(pvalueAdj = p.adjust(pvalue, method = adjust)) %>%
      column_to_rownames("ID")
    
    non_param_mann <- bind_cols(group_means, non_param_mann) %>%
      rownames_to_column("ID") %>%
      mutate(Fold_Change_Ratio = as.numeric(round(group_means[,2]/group_means[,1], 3)),
             Difference_Of_Means = as.numeric(round(group_means[,2] - group_means[,1], 3))) %>% 
      column_to_rownames("ID") %>%
      select(1,2,5,6,3,4)

    return(non_param_mann)
  }

  else if (method == "kruskal"){

    non_param_kru <- data.frame(pvalue = apply(e, 2, function(x){kruskal.test(x ~ as.factor(Group))$p.value}))
    non_param_kru <- non_param_kru %>%
      mutate(pvalueAdj = p.adjust(pvalue, method = adjust),
             Kruskal_Wallis_Rank_Sum = apply(e, 2, function(x){kruskal.test(x ~ as.factor(Group))$statistic}))

    non_param_kru <- bind_cols(group_means, non_param_kru)
    
    return(non_param_kru)
  }

}

