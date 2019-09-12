
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

    model <- lmFit(trans_limma, initialmodel)
    model <- contrasts.fit(model, cont.matrix)

    modelstats <- eBayes(model)
    res <- topTable(modelstats, number = ncol(data_limma),
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

    model2 <- lmFit(trans_limma2, initialmodel2)
    model2 <- contrasts.fit(model2, cont.matrix2)

    modelstats2 <- eBayes(model2)
    res2 <- topTable(modelstats2, number= ncol(data_limma) ,
                     coef = contrast, sort.by = "p", adjust.method = adjust)

    return(res2)
  }

}

