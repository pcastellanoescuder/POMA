
#' Implementation of Limma R Package in Metabolomics
#'
#' @description PomaLimma() uses the classical limma package for metabolomics.
#'
#' @param data_limma A data frame with metabolites. First column must be the subject ID and second column must be a factor with the subject group.
#' @param contrast A character with the limma comparison. For example, "Group1-Group2" or "control-intervention".
#' @param covariates A data frame with covariates. The first column must be the subject ID in the same order as in the metabolites data (optional).
#' @param adjust Multiple comparisons correction method.
#'
#' @export
#'
#' @return A data frame with the limma results.
#' @references Matthew E. Ritchie, Belinda Phipson, Di Wu, Yifang Hu, Charity W. Law, Wei Shi, Gordon K. Smyth, limma powers differential expression analyses for RNA-sequencing and microarray studies, Nucleic Acids Research, Volume 43, Issue 7, 20 April 2015, Page e47, https://doi.org/10.1093/nar/gkv007
#' @author Pol Castellano-Escuder
PomaLimma <- function(data_limma,
                      contrast = NULL,
                      covariates = NULL,
                      adjust = c("fdr", "holm", "hochberg", "hommel", "bonferroni", "BH", "BY")){

  if (is.null(contrast)) {
    stop(crayon::red(clisymbols::symbol$cross, "Contrast argument is empty! You have to specify a contrast."))
  }
  if (missing(adjust)) {
    adjust <- "fdr"
    warning("adjust argument is empty! FDR will be used")
  }

  colnames(data_limma)[1:2] <- c("ID", "Group")

  contrasts <- levels(as.factor(data_limma$Group))
  fac1 <- as.factor(data_limma$Group)

  if (is.null(covariates)){

    initialmodel <- model.matrix( ~ 0 + fac1)
    colnames(initialmodel) <- contrasts

    cont.matrix <- limma::makeContrasts(contrasts = contrast,
                                        levels = initialmodel)

    trans_limma <- t(data_limma[, c(3:ncol(data_limma))])

    model <- limma::lmFit(trans_limma, initialmodel)
    model <- limma::contrasts.fit(model, cont.matrix)

    modelstats <- limma::eBayes(model)
    res <- limma::topTable(modelstats, number = ncol(data_limma),
                           coef = contrast, sort.by = "p", adjust.method = adjust)

    return(res)

  }

  ####

  else if (!is.null(covariates)){

    colnames(covariates)[1] <- "ID"

    form <- as.formula(noquote(paste("~ 0 + fac1 + ",
                                     paste0(colnames(covariates)[2:length(covariates)],
                                            collapse = " + ", sep=""), sep = "")))

    initialmodel2 <- model.matrix(form , covariates)
    colnames(initialmodel2)[1:length(levels(fac1))] <- contrasts

    cont.matrix2 <- limma::makeContrasts(contrasts = contrast,
                                         levels = initialmodel2)

    trans_limma2 <- t(data_limma[, c(3:ncol(data_limma))])

    model2 <- limma::lmFit(trans_limma2, initialmodel2)
    model2 <- limma::contrasts.fit(model2, cont.matrix2)

    modelstats2 <- limma::eBayes(model2)
    res2 <- limma::topTable(modelstats2, number= ncol(data_limma) ,
                            coef = contrast, sort.by = "p", adjust.method = adjust)

    return(res2)
  }

}

