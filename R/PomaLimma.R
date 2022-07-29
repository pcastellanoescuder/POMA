
#' Easy Implementation of `limma` Bioconductor Package
#'
#' @description PomaLimma() uses the classical `limma` package.
#'
#' @param data A SummarizedExperiment object.
#' @param contrast A character with the limma comparison. For example, "Group1-Group2" or "control-intervention".
#' @param covariates Logical. If it's set to `TRUE` all metadata variables stored in `colData` will be used as covariables. Default = FALSE.
#' @param covs Character vector indicating the name of `colData` columns that will be included as covariates. Default is NULL (all variables).
#' @param adjust Multiple comparisons correction method. Options are: "fdr", "holm", "hochberg", "hommel", "bonferroni", "BH" and "BY".
#' @param cutoff Default is NULL. If this value is replaced for a numeric value, the resultant table will contains only those features with an adjusted p-value below selected value. 
#'
#' @export
#'
#' @return A tibble with limma results.
#' @references Matthew E. Ritchie, Belinda Phipson, Di Wu, Yifang Hu, Charity W. Law, Wei Shi, Gordon K. Smyth, limma powers differential expression analyses for RNA-sequencing and microarray studies, Nucleic Acids Research, Volume 43, Issue 7, 20 April 2015, Page e47, https://doi.org/10.1093/nar/gkv007
#' @author Pol Castellano-Escuder
#'
#' @importFrom magrittr %>%
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
                      covs = NULL,
                      adjust = "fdr",
                      cutoff = NULL){

  if (missing(data)) {
    stop("data argument is empty!")
  }
  if(!is(data[1], "SummarizedExperiment")){
    stop("data is not a SummarizedExperiment object. \nSee POMA::PomaSummarizedExperiment or SummarizedExperiment::SummarizedExperiment")
  }
  if (is.null(contrast)) {
    stop("Contrast argument is empty! Specify a contrast.")
  }
  if (!(adjust %in% c("fdr", "holm", "hochberg", "hommel", "bonferroni", "BH", "BY"))) {
    stop("Incorrect value for adjust argument!")
  }
  if(covariates & ncol(SummarizedExperiment::colData(data)) == 1){
    stop("Seems there aren't covariates in your data...")
  }

  Group <- SummarizedExperiment::colData(data)[,1]
  e <- SummarizedExperiment::assay(data)

  contrasts <- levels(as.factor(Group))
  fac1 <- as.factor(Group)

  if (!covariates){

    initialmodel <- stats::model.matrix( ~ 0 + fac1)
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
    
    res <- res %>% 
      tibble::rownames_to_column("feature") %>% 
      dplyr::as_tibble()
      
    return(res)

  }

  ##

  else {

    if(is.null(covs)){
      covariates <- SummarizedExperiment::colData(data) %>%
        as.data.frame() %>%
        dplyr::select(-1)
    } 
    else{
      covariates <- SummarizedExperiment::colData(data) %>%
        as.data.frame() %>%
        dplyr::select(-1) %>% 
        dplyr::select_at(dplyr::vars(dplyr::matches(covs)))
    }
    
    form <- as.formula(noquote(paste("~ 0 + fac1 + ", paste0(colnames(covariates), collapse = " + ", sep=""), sep = "")))

    initialmodel2 <- stats::model.matrix(form, covariates)
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
    
    res2 <- res2 %>% 
      tibble::rownames_to_column("feature") %>% 
      dplyr::as_tibble()
    
    return(res2)
    
  }

}

