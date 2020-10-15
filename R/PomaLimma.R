
#' Implementation of limma R Package on Mass Spectrometry Data
#'
#' @description PomaLimma() uses the classical limma package for MS data.
#'
#' @param data A MSnSet object. First `pData` column must be the subject group/type.
#' @param contrast A character with the limma comparison. For example, "Group1-Group2" or "control-intervention".
#' @param covariates Logical. If it's set to `TRUE` all metadata variables stored in `pData` will be used as covariables. Default = FALSE.
#' @param adjust Multiple comparisons correction method. Options are: "fdr", "holm", "hochberg", "hommel", "bonferroni", "BH" and "BY".
#' @param cutoff Default is NULL. If this value is replaced for a numeric value, the resultant table will contains only those features with an adjusted p-value below selected value. 
#'
#' @export
#'
#' @return A data frame with limma results.
#' @references Matthew E. Ritchie, Belinda Phipson, Di Wu, Yifang Hu, Charity W. Law, Wei Shi, Gordon K. Smyth, limma powers differential expression analyses for RNA-sequencing and microarray studies, Nucleic Acids Research, Volume 43, Issue 7, 20 April 2015, Page e47, https://doi.org/10.1093/nar/gkv007
#' @author Pol Castellano-Escuder
#'
#' @importFrom limma makeContrasts lmFit contrasts.fit eBayes topTable
#' @importFrom crayon red
#' @importFrom magrittr %>%
#' @importFrom dplyr filter
#' @importFrom clisymbols symbol
#' @importFrom Biobase varLabels pData exprs
#' 
#' @examples 
#' data("st000284")
#' 
#' st000284 %>%
#'   PomaNorm() %>%
#'   PomaLimma(contrast = "Healthy-CRC", adjust = "fdr")
PomaLimma <- function(data,
                      contrast = NULL,
                      covariates = FALSE,
                      adjust = "fdr",
                      cutoff = NULL){

  if (missing(data)) {
    stop(crayon::red(clisymbols::symbol$cross, "data argument is empty!"))
  }
  if(!is(data[1], "MSnSet")){
    stop(paste0(crayon::red(clisymbols::symbol$cross, "data is not a MSnSet object."), 
                " \nSee POMA::PomaMSnSetClass or MSnbase::MSnSet"))
  }
  if (is.null(contrast)) {
    stop(crayon::red(clisymbols::symbol$cross, "Contrast argument is empty! You have to specify a contrast."))
  }
  if (!(adjust %in% c("fdr", "holm", "hochberg", "hommel", "bonferroni", "BH", "BY"))) {
    stop(crayon::red(clisymbols::symbol$cross, "Incorrect value for adjust argument!"))
  }
  if (missing(adjust)) {
    warning("adjust argument is empty! FDR will be used")
  }
  if(isTRUE(covariates) & ncol(pData(data)) == 1){
    stop(crayon::red(clisymbols::symbol$cross, "Seems that your data don't have covariates..."))
  }

  Biobase::varLabels(data)[1] <- "Group"
  Group <- Biobase::pData(data)$Group
  e <- Biobase::exprs(data)

  contrasts <- levels(as.factor(Group))
  fac1 <- as.factor(Group)

  if (!covariates){

    initialmodel <- model.matrix( ~ 0 + fac1)
    colnames(initialmodel) <- contrasts

    cont.matrix <- limma::makeContrasts(contrasts = contrast,
                                        levels = initialmodel)

    model <- limma::lmFit(e, initialmodel)
    model <- limma::contrasts.fit(model, cont.matrix)

    modelstats <- limma::eBayes(model)
    res <- limma::topTable(modelstats, number = nrow(e),
                           coef = contrast, sort.by = "p", adjust.method = adjust)

    if(!is.null(cutoff)){
      res <- res %>%
        dplyr::filter(adj.P.Val <= cutoff)
    }
    
    return(res)

  }

  ##

  else {

    covariates <- as.data.frame(pData(data)[, 2:ncol(pData(data))])
    
    if(ncol(covariates) == 1){
      colnames(covariates) <- colnames(pData(data))[2]
    }

    form <- as.formula(noquote(paste("~ 0 + fac1 + ", paste0(colnames(covariates), collapse = " + ", sep=""), sep = "")))

    initialmodel2 <- model.matrix(form, covariates)
    colnames(initialmodel2)[1:length(levels(fac1))] <- contrasts

    cont.matrix2 <- limma::makeContrasts(contrasts = contrast,
                                         levels = initialmodel2)

    model2 <- limma::lmFit(e, initialmodel2)
    model2 <- limma::contrasts.fit(model2, cont.matrix2)

    modelstats2 <- limma::eBayes(model2)
    res2 <- limma::topTable(modelstats2, number = nrow(e) ,
                            coef = contrast, sort.by = "p", adjust.method = adjust)

    if(!is.null(cutoff)){
      res2 <- res2 %>%
        dplyr::filter(adj.P.Val <= cutoff)
    }
    
    return(res2)
    
  }

}

