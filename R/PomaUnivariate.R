
#' Univariate Statistical Methods for Mass Spectrometry Data
#'
#' @description PomaUnivariate() allows users to perform different univariate statistical analysis on MS data.
#'
#' @param data_uni A MSnSet object. First `pData` column must be the suject group/type.
#' @param covariates Logical. If it's set to `TRUE` all metadata variables stored in `pData` will be used as covariables. Default = FALSE.
#' @param method Univariate statistical method. Options are c("ttest", "anova", "mann", "kruskal").
#' @param paired Logical indicates if the data is paired or not.
#' @param var_equal Logical indicates if the data variance is equal or not.
#' @param adjust Multiple comparisons correction method.
#'
#' @export
#'
#' @return A data frame with results.
#' @author Pol Castellano-Escuder
#'
#' @importFrom tibble rownames_to_column column_to_rownames as_tibble
#' @importFrom dplyr select mutate
#' @importFrom magrittr %>%
#' @importFrom crayon red
#' @importFrom clisymbols symbol
#' @importFrom Biobase varLabels pData exprs
PomaUnivariate <- function(data_uni,
                           covariates = FALSE,
                           method = c("ttest", "anova", "mann", "kruskal"),
                           paired = FALSE,
                           var_equal = FALSE,
                           adjust = c("fdr", "holm", "hochberg", "hommel", "bonferroni", "BH", "BY")){

  if (missing(method)) {
    stop(crayon::red(clisymbols::symbol$cross, "Select a method!"))
  }
  if (!(method %in% c("ttest", "anova", "mann", "kruskal"))) {
    stop(crayon::red(clisymbols::symbol$cross, "Incorrect value for method argument!"))
  }
  if (missing(adjust)) {
    adjust <- "fdr"
    warning("adjust argument is empty! FDR will be used")
  }

  Biobase::varLabels(data_uni)[1] <- "Group"
  Group <- Biobase::pData(data_uni)$Group
  e <- t(Biobase::exprs(data_uni))

  if(method == "ttest"){

    stat <- function(x){t.test(x ~ Group, na.rm = TRUE, alternative = c("two.sided"),
                               var.equal = eval(parse(text = var_equal)),
                               paired = eval(parse(text = paired)))$p.value}
    stat_G2 <- function(x){t.test(x ~ Group, na.rm = TRUE, alternative = c("two.sided"),
                                  var.equal = eval(parse(text = var_equal)))$estimate[[2]]}
    stat_G1 <- function(x){t.test(x ~ Group, na.rm = TRUE, alternative = c("two.sided"),
                                  var.equal = eval(parse(text = var_equal)))$estimate[[1]]}

    p <- data.frame(pvalue = apply(FUN = stat, MARGIN = 2, X = e))

    p <- p %>%
      rownames_to_column("metabolite") %>%
      as_tibble() %>%
      mutate(pvalue_Adj = p.adjust(pvalue, method = adjust)) %>%
      column_to_rownames("metabolite")

    G2 <- round(data.frame(Mean_G2 = apply(FUN = stat_G2, MARGIN = 2, X = e)), 3)
    G1 <- round(data.frame(Mean_G1 = apply(FUN = stat_G1, MARGIN = 2, X = e)), 3)
    means <- cbind(G1, G2)

    means <- means %>%
      rownames_to_column("metabolite") %>%
      mutate(Fold_Change_Ratio = as.numeric(round(Mean_G2/Mean_G1, 3)),
             Difference_Of_Means = as.numeric(round(Mean_G1 - Mean_G2, 3))) %>%
      column_to_rownames("metabolite")

    p <- cbind(means, p)

    return(p)
  }

  else if(method == "anova"){

    if(!isTRUE(covariates)){

      stat2 <- function(x){anova(aov(x ~ Group))$"Pr(>F)"[1]}
      p2 <- data.frame(pvalue = apply(FUN = stat2, MARGIN = 2, X = e))

      p2 <- p2 %>%
        rownames_to_column("metabolite") %>%
        as_tibble() %>%
        mutate(pvalue_Adj = p.adjust(pvalue, method = adjust)) %>%
        column_to_rownames("metabolite")

      return(p2)

    }
    else{

      covariate_uni <- merge(e, pData(data_uni), by = "row.names")
      covariate_uni <- covariate_uni %>% select(-Row.names, -Group)

      stat3 <- function(x){anova(aov(x ~ Group))$"Pr(>F)"[1]}
      p3 <- data.frame(apply(FUN = stat3, MARGIN = 2, X = covariate_uni))

      colnames(p3) <- "pvalue"

      p3 <- p3 %>%
        rownames_to_column %>%
        mutate(pvalue_Adj = p.adjust(pvalue, method = adjust)) %>%
        column_to_rownames("rowname")

      return(p3)

    }
  }

  else if(method == "mann"){

    non_param_mann <- data.frame(pvalue = apply(e, 2,
                                                 function(x){wilcox.test(x ~ as.factor(Group),
                                                                         paired = eval(parse(text = paired)))$p.value}))

    non_param_mann <- non_param_mann %>%
      rownames_to_column %>%
      mutate(pvalue_Adj = p.adjust(pvalue, method = adjust)) %>%
      column_to_rownames("rowname")

    manndata <- data.frame(cbind(Group = Group, e))
    manndata[,2:ncol(manndata)] <- sapply(manndata[,2:ncol(manndata)], as.numeric)

    means <- aggregate(manndata[, 2:ncol(manndata)], list(manndata$Group), mean)

    means <- column_to_rownames(means, var = "Group.1")
    means <- as_tibble(t(means), rownames = NA)
    colnames(means) <- c("Mean_G1", "Mean_G2")

    means <- means %>%
      rownames_to_column("metabolite") %>%
      mutate(Fold_Change_Ratio = as.numeric(round(Mean_G2/Mean_G1, 3)),
             Difference_Of_Means = as.numeric(round(Mean_G1 - Mean_G2, 3))) %>%
      column_to_rownames("metabolite")

    non_param_mann <- cbind(means, non_param_mann)

    return(non_param_mann)
  }

  else if (method == "kruskal"){

    non_param_kru <- data.frame(pvalue = apply(e, 2, function(x){kruskal.test(x ~ as.factor(Group))$p.value}))
    non_param_kru <- non_param_kru %>%
      rownames_to_column("metabolite") %>%
      as_tibble() %>%
      mutate(pvalue_Adj = p.adjust(pvalue, method = adjust),
             Kruskal_Wallis_Rank_Sum = apply(e, 2, function(x){kruskal.test(x ~ as.factor(Group))$statistic})) %>%
      column_to_rownames("metabolite")

    return(non_param_kru)
  }

}

