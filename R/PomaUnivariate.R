
#' Univariate Statistical Methods for Metabolomics
#'
#' @description PomaUnivariate() allows users to perform different univariate statistical analysis on metabolomic data.
#'
#' @param data_uni A data frame with metabolites. First column must be the subject ID and second column must be a factor with the subject group.
#' @param covariates A data frame with covariates. The first column must be the subject ID in the same order as in the metabolites data (optional).
#' @param method Univariate statistical method. Options are c("ttest", "anova", "mann", "kruskal").
#' @param paired Logical indicates if the data is paired or not.
#' @param var_equal Logical indicates if the data variance is equal or not.
#' @param adjust Multiple comparisons correction method.
#'
#' @export
#'
#' @return A data frame with the results.
#' @author Pol Castellano-Escuder
PomaUnivariate <- function(data_uni,
                           covariates = NULL,
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


  colnames(data_uni)[1:2] <- c("ID", "Group")
  Group <- data_uni$Group

  if(method == "ttest"){

    stat <- function(x){t.test(x ~ Group, na.rm = TRUE, alternative = c("two.sided"),
                               var.equal = eval(parse(text = var_equal)),
                               paired = eval(parse(text = paired)))$p.value}
    stat_G2 <- function(x){t.test(x ~ Group, na.rm = TRUE, alternative = c("two.sided"),
                                  var.equal = eval(parse(text = var_equal)))$estimate[[2]]}
    stat_G1 <- function(x){t.test(x ~ Group, na.rm = TRUE, alternative = c("two.sided"),
                                  var.equal = eval(parse(text = var_equal)))$estimate[[1]]}

    p <- data.frame(P.Value = apply(FUN = stat, MARGIN = 2, X = data_uni[,c(3:ncol(data_uni))]))

    p <- p %>%
      rownames_to_column("metabolite") %>%
      as_tibble() %>%
      mutate(adj.P.Val = p.adjust(P.Value, method = adjust)) %>%
      column_to_rownames("metabolite")

    G2 <- round(data.frame(Mean_G2 = apply(FUN = stat_G2, MARGIN = 2, X = data_uni[,c(3:ncol(data_uni))])), 3)
    G1 <- round(data.frame(Mean_G1 = apply(FUN = stat_G1, MARGIN = 2, X = data_uni[,c(3:ncol(data_uni))])), 3)
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

    if(is.null(covariates)){

      stat2 <- function(x){anova(aov(x ~ Group, data = data_uni))$"Pr(>F)"[1]}
      p2 <- data.frame(P.Value = apply(FUN = stat2, MARGIN = 2, X = data_uni[,c(3:ncol(data_uni))]))

      p2 <- p2 %>%
        rownames_to_column("metabolite") %>%
        as_tibble() %>%
        mutate(adj.P.Val = p.adjust(P.Value, method = adjust)) %>%
        column_to_rownames("metabolite")

      return(p2)

    }
    else{

      colnames(covariates)[1] <- "ID"
      covariate_uni <- merge(data_uni, covariates, by = "ID")

      form2 <- noquote(paste("y ~", paste(colnames(covariate_uni)[c(2,(length(data_uni)+1):length(covariate_uni))],
                                          collapse = " + ", sep="")))

      stat3 <- function(y){anova(aov(as.formula(form2), data = covariate_uni))$"Pr(>F)"}
      p3 <- data.frame(apply(FUN = stat3, MARGIN = 2, X = covariate_uni[,3:length(data_uni)]))

      p3 <- p3 %>%
        rownames_to_column %>%
        as_tibble() %>%
        gather(var, value, -rowname) %>%
        spread(rowname, value) %>%
        column_to_rownames("var") %>%
        setNames(paste0("P.Value_", colnames(covariates)))

      p3 <- p3[, 1:ncol(p3)-1] %>%
        rename(P.Value_Group = "P.Value_ID")

      p3 <- p3 %>%
        rownames_to_column("metabolite") %>%
        as_tibble() %>%
        mutate_at(vars(contains("P.Value")), list(adjusted = ~ p.adjust(. , method = adjust))) %>%
        column_to_rownames("metabolite")

      return(p3)

    }
  }

  else if(method == "mann"){

    non_param_mann <- data.frame(P.Value = apply(data_uni[,3:ncol(data_uni)], 2,
                                                 function(x){wilcox.test(x ~ as.factor(Group),
                                                                         paired = eval(parse(text = paired)))$p.value}))

    non_param_mann <- non_param_mann %>% as_tibble() %>% mutate(adj.P.Val = p.adjust(P.Value, method = adjust))

    means <- data_uni %>% group_by(Group) %>% summarise_at(vars(names(data_uni[3:ncol(data_uni)])), mean)
    means <- column_to_rownames(means, var = "Group")
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

    non_param_kru <- data.frame(P.Value = apply(data_uni[,3:ncol(data_uni)], 2, function(x){kruskal.test(x ~ as.factor(Group))$p.value}))
    non_param_kru <- non_param_kru %>%
      rownames_to_column("metabolite") %>%
      as_tibble() %>%
      mutate(adj.P.Val = p.adjust(P.Value, method = adjust),
             Kruskal_Wallis_Rank_Sum = apply(data_uni[,3:ncol(data_uni)], 2, function(x){kruskal.test(x ~ as.factor(Group))$statistic})) %>%
      column_to_rownames("metabolite")

    return(non_param_kru)
  }

}

