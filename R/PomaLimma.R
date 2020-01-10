
#' Implementation of Limma R Package in Mass Spectrometry Data
#'
#' @description PomaLimma() uses the classical limma package for MS data.
#'
#' @param data_limma A MSnSet object. First `pData` column must be the suject group/type.
#' @param contrast A character with the limma comparison. For example, "Group1-Group2" or "control-intervention".
#' @param covariates Logical. If it's set to `TRUE` all metadata variables stored in `pData` will be used as covariables. Default = FALSE.
#' @param adjust Multiple comparisons correction method.
#'
#' @export
#'
#' @return A data frame with limma results.
#' @references Matthew E. Ritchie, Belinda Phipson, Di Wu, Yifang Hu, Charity W. Law, Wei Shi, Gordon K. Smyth, limma powers differential expression analyses for RNA-sequencing and microarray studies, Nucleic Acids Research, Volume 43, Issue 7, 20 April 2015, Page e47, https://doi.org/10.1093/nar/gkv007
#' @author Pol Castellano-Escuder
#'
#' @importFrom limma makeContrasts lmFit contrasts.fit eBayes topTable
#' @importFrom crayon red
#' @importFrom clisymbols symbol
#' @importFrom Biobase varLabels pData exprs
PomaLimma <- function(data_limma,
                      contrast = NULL,
                      covariates = FALSE,
                      adjust = c("fdr", "holm", "hochberg", "hommel", "bonferroni", "BH", "BY")){

  if (is.null(contrast)) {
    stop(crayon::red(clisymbols::symbol$cross, "Contrast argument is empty! You have to specify a contrast."))
  }
  if (missing(adjust)) {
    adjust <- "fdr"
    warning("adjust argument is empty! FDR will be used")
  }

  Biobase::varLabels(data_limma)[1] <- "Group"
  Group <- Biobase::pData(data_limma)$Group
  e <- Biobase::exprs(data_limma)

  contrasts <- levels(as.factor(Group))
  fac1 <- as.factor(Group)

  if (isFALSE(covariates)){

    initialmodel <- model.matrix( ~ 0 + fac1)
    colnames(initialmodel) <- contrasts

    cont.matrix <- limma::makeContrasts(contrasts = contrast,
                                        levels = initialmodel)

    model <- limma::lmFit(e, initialmodel)
    model <- limma::contrasts.fit(model, cont.matrix)

    modelstats <- limma::eBayes(model)
    res <- limma::topTable(modelstats, number = nrow(e),
                           coef = contrast, sort.by = "p", adjust.method = adjust)

    return(res)

  }

  ####

  else {

    covariates <- pData(data_limma)[, 2:ncol(pData(data_limma))]

    form <- as.formula(noquote(paste("~ 0 + fac1 + ", paste0(colnames(covariates), collapse = " + ", sep=""), sep = "")))

    initialmodel2 <- model.matrix(form , covariates)
    colnames(initialmodel2)[1:length(levels(fac1))] <- contrasts

    cont.matrix2 <- limma::makeContrasts(contrasts = contrast,
                                         levels = initialmodel2)

    model2 <- limma::lmFit(e, initialmodel2)
    model2 <- limma::contrasts.fit(model2, cont.matrix2)

    modelstats2 <- limma::eBayes(model2)
    res2 <- limma::topTable(modelstats2, number = nrow(e) ,
                            coef = contrast, sort.by = "p", adjust.method = adjust)

    return(res2)
  }

}

