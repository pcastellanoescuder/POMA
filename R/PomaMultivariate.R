
#' Multivariate Statistical Methods for Mass Spectrometry Data
#'
#' @description PomaMultivariate() allows users to perform different multivariate statistical analysis on MS data.
#'
#' @param data_multi A MSnSet object. First `pData` column must be the suject group/type.
#' @param method A multivariate method. Options are c("pca", "plsda", "splsda").
#' @param components Numeric. Number of components to include in the model. Default is 5.
#' @param center Logical that indicates whether the variables should be shifted to be zero centered. Default is FALSE.
#' @param scale Logical that indicates whether the variables should be scaled to have unit variance before the analysis takes place. Default is FALSE.
#' @param ellipse Logical that indicates whether a 95%CI ellipse should be plotted in scores plot. Default is TRUE.
#' @param validation (Only for "plsda" and "splsda" methods) Validation method. Options are c("Mfold", "loo").
#' @param folds (Only for "plsda" and "splsda" methods) Numeric. Number of folds for Mfold validation method (default is 5). If the validation method is loo, this value will become to 1.
#' @param nrepeat (Only for "plsda" and "splsda" methods) Numeric. Number of iterations for the validation method selected.
#' @param num_features (Only for "splsda" method) Numeric. Number of variables selected to discriminate groups.
#'
#' @export
#'
#' @return A list with all results for multivariate statistical analysis including plots and data frames.
#' @author Pol Castellano-Escuder
#'
#' @import ggplot2
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr select
#' @importFrom magrittr %>%
#' @importFrom mixOmics pca plsda perf vip tune.splsda splsda selectVar
#' @importFrom reshape2 melt
#' @importFrom crayon red
#' @importFrom clisymbols symbol
#' @importFrom Biobase varLabels pData exprs
PomaMultivariate <- function(data_multi,
                             method = c("pca", "plsda", "splsda"),
                             components = 5,
                             center = FALSE,
                             scale = FALSE,
                             ellipse = TRUE,
                             validation = c("Mfold", "loo"),
                             folds = 5,
                             nrepeat = 10,
                             num_features = 10){

  if (missing(method)) {
    stop(crayon::red(clisymbols::symbol$cross, "Select a method!"))
  }
  if (!(method %in% c("pca", "plsda", "splsda"))) {
    stop(crayon::red(clisymbols::symbol$cross, "Incorrect value for method argument!"))
  }
  if (missing(validation)) {
    validation <- "Mfold"
    if (method %in% c("plsda", "splsda")){
      warning("validation argument is empty! Mfold will be used")
    }
  }
  if (!(validation %in% c("Mfold", "loo"))) {
    stop(crayon::red(clisymbols::symbol$cross, "Incorrect validation method! Please choose 'Mfold' or 'loo'"))
  }

  Biobase::varLabels(data_multi)[1] <- "Group"

  df <- t(Biobase::exprs(data_multi))

  if(method == "pca"){

    X <- as.matrix(df)
    Y <- as.factor(Biobase::pData(data_multi)$Group)
    pca_res <- mixOmics::pca(X, ncomp = components, center = center, scale = scale)

    PCi <- data.frame(pca_res$x, Groups = Y)

    scoresplot <- ggplot(PCi, aes(x = PC1, y = PC2, col = Groups))+
      geom_point(size = 3, alpha = 0.5) +
      xlab(paste0("PC1 (", round(100*(pca_res$explained_variance)[1], 2), "%)")) +
      ylab(paste0("PC2 (", round(100*(pca_res$explained_variance)[2], 2), "%)")) +
      {if(ellipse)stat_ellipse(type = "norm")} +
      theme_bw()

    ####

    eigenvalues <- data.frame(Percent_Variance_Explained = round(pca_res$explained_variance*100, 4))
    eigenvalues <- eigenvalues %>% rownames_to_column("Principal_Component")

    screeplot <- ggplot(eigenvalues, aes(x = Principal_Component, y = Percent_Variance_Explained, fill = NULL)) +
      geom_bar(stat = "identity", fill = rep(c("lightblue"), nrow(eigenvalues))) +
      xlab("Principal Component") +
      ylab("% Variance Explained") +
      geom_label(data = eigenvalues, aes(label =  paste0(round(Percent_Variance_Explained, 3), "%"))) +
      theme_bw()

    ####

    score_data <- PCi %>% dplyr::select(-Groups) %>% round(4)

    ####

    pca_res2 <- mixOmics::pca(X, ncomp = components, center = T, scale = T)

    PCi2 <- data.frame(pca_res2$x, Groups = Y)

    PCAloadings <- data.frame(pca_res2$loadings$X)

    biplot <- ggplot(PCi2, aes(x = PC1, y = PC2, col = Groups))+
      geom_point(size = 3, alpha = 0.5) +
      xlab(paste0("PC1 (", round(100*(pca_res2$explained_variance)[1], 2), "%)")) +
      ylab(paste0("PC2 (", round(100*(pca_res2$explained_variance)[2], 2), "%)")) +
      theme_bw() +
      {if(ellipse)stat_ellipse(type = "norm")} +
      geom_segment(data = PCAloadings,
                   aes(x = 0, y = 0,
                       xend = (PC1*12),
                       yend = (PC2*12)),
                   arrow = arrow(length = unit(1/2, "picas")), color = "grey19") +
      annotate("text", x = (PCAloadings$PC1*12),
               y = (PCAloadings$PC2*12),
               label = rownames(PCAloadings), size = 4)

    return(list(screeplot = screeplot, scoresplot = scoresplot,
                score_data = score_data, eigenvalues = eigenvalues, biplot = biplot))

  }

  else if (method == "plsda"){

    X <- as.matrix(df)
    Y <- as.factor(Biobase::pData(data_multi)$Group)

    plsda_res <- mixOmics::plsda(X, Y, ncomp = components)

    PLSDAi <- data.frame(plsda_res$variates$X, Groups = Y)

    scoresplot <- ggplot(PLSDAi, aes(x = comp1, y = comp2, col = Groups))+
      geom_point(size=3,alpha=0.5) +
      xlab("Component 1") +
      ylab("Component 2") +
      {if(ellipse)stat_ellipse(type = "norm")} +
      theme_bw()

    #####

    set.seed(69)
    perf_plsda <- mixOmics::perf(plsda_res, validation = validation, folds = folds,
                                progressBar = TRUE, auc = TRUE, nrepeat = nrepeat)

    overall <- data.frame(perf_plsda$error.rate[1]) %>%
      round(4) %>%
      rownames_to_column("Component")

    ber <- data.frame(perf_plsda$error.rate[2]) %>%
      round(4) %>%
      rownames_to_column("Component")

    errors_plsda1 <- reshape2::melt(ber, id.vars=c("Component"))
    errors_plsda2 <- reshape2::melt(overall, id.vars=c("Component"))
    errors_plsda <- rbind(errors_plsda1, errors_plsda2)

    errors_plsda_plot <- ggplot(data = errors_plsda, aes(x = Component, y = value, group = variable)) +
      geom_line(aes(color=variable)) +
      geom_point(aes(color=variable)) +
      ylab("Error") +
      theme_bw()

    ####

    plsda_vip <- data.frame(mixOmics::vip(plsda_res))
    plsda_vip <- plsda_vip[order(plsda_vip[,1], decreasing = T) ,]

    plsda_vip <- plsda_vip %>% rownames_to_column("Variable")
    plsda_vip_top <- plsda_vip[1:15 ,]

    plsda_vip_top <- plsda_vip_top %>%
      mutate(Variable = factor(Variable, levels = Variable[order(comp1)]))

    vip_plsda_plot <- ggplot(plsda_vip_top, aes(x = Variable, y = comp1, fill = NULL)) +
      geom_bar(stat="identity", fill = rep(c("lightblue"), nrow(plsda_vip_top))) +
      coord_flip() +
      ylab("VIP") +
      xlab("") +
      geom_label(data = plsda_vip_top, aes(label = round(comp1, 2))) +
      theme_bw()

    ####

    scores_plsda <- PLSDAi %>% dplyr::select(-Groups) %>% round(4)

    return(list(scoresplot = scoresplot, errors_plsda = errors_plsda,
                errors_plsda_plot = errors_plsda_plot, plsda_vip_table = plsda_vip,
                vip_plsda_plot = vip_plsda_plot, scores_plsda = scores_plsda))
  }

  else if (method == "splsda"){

    X <- as.matrix(df)
    Y <- as.factor(Biobase::pData(data_multi)$Group)

    list_keepX <- c(1:num_features)

    tune_splsda <- mixOmics::tune.splsda(X, Y, ncomp = components, validation = validation, folds = folds,
                                         progressBar = TRUE, dist = 'max.dist', measure = "BER",
                                         test.keepX = list_keepX, nrepeat = nrepeat) # cpus = 4

    error <- tune_splsda$error.rate
    ncomp <- tune_splsda$choice.ncomp$ncomp # optimal number of components based on t-tests
    select_keepX <- tune_splsda$choice.keepX[1:ncomp]  # optimal number of variables to select

    errors_splsda_out <- data.frame(tune_splsda$error.rate) %>% round(4) %>%
      rownames_to_column("features")
    errors_splsda <- reshape2::melt(errors_splsda_out, id.vars=c("features"))

    errors_sd <- data.frame(tune_splsda$error.rate.sd) %>% rownames_to_column("features_sd")
    errors_sd <- reshape2::melt(errors_sd, id.vars=c("features_sd"))

    errors_splsda <- cbind(errors_splsda, sd = errors_sd[,3])

    bal_error_rate <- ggplot(data = errors_splsda, aes(x = features, y = value, group = variable)) +
      geom_line(aes(color=variable)) +
      geom_point(aes(color=variable)) +
      geom_errorbar(aes(ymin = value-sd, ymax = value+sd), width=.1) +
      theme_bw() +
      xlab("Number of features") +
      ylab("Error")

    ####

    if (ncomp == 1){
      ncompX <- 2
    }else{
      ncompX <- ncomp}

    ####

    res_splsda <- mixOmics::splsda(X, Y, ncomp = ncompX, keepX = select_keepX)

    SPLSDAi <- data.frame(res_splsda$variates$X, Groups = Y)

    splsda_scores_plot <- ggplot(SPLSDAi, aes(x = comp1, y = comp2, col = Groups)) +
      geom_point(size=3,alpha=0.5) +
      xlab("Component 1") +
      ylab("Component 2") +
      {if(ellipse)stat_ellipse(type = "norm")} +
      theme_bw()

    scores_splsda <- SPLSDAi %>% dplyr::select(-Groups) %>% round(4)

    selected_variables <- mixOmics::selectVar(res_splsda, comp = 1)
    selected_variables <- round(selected_variables$value, 4)
    selected_variables <- data.frame(Feature = rownames(selected_variables), Value = selected_variables$value.var)

    ####

    return(list(ncomp = ncomp, select_keepX = select_keepX, errors_splsda = errors_splsda,
                bal_error_rate = bal_error_rate, splsda_scores_plot = splsda_scores_plot,
                scores_splsda = scores_splsda, selected_variables = selected_variables))

  }

}

