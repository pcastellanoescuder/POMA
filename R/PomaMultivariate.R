
#' Multivariate Statistical Methods for Mass Spectrometry Data
#'
#' @description PomaMultivariate() allows users to perform different multivariate statistical analysis on MS data.
#'
#' @param data A MSnSet object. First `pData` column must be the subject group/type.
#' @param method A multivariate method. Options are: "pca", "plsda" and "splsda".
#' @param components Numeric. Number of components to include in the model. Default is 5.
#' @param center Logical that indicates whether the variables should be shifted to be zero centered. Default is FALSE.
#' @param scale Logical that indicates whether the variables should be scaled to have unit variance before the analysis takes place. Default is FALSE.
#' @param labels Logical indicating if sample names should be plotted or not.
#' @param load_length Numeric between 1 and 2. Define the length of biplot loadings. Default is 1.
#' @param ellipse Logical that indicates whether a 95%CI ellipse should be plotted in scores plot. Default is TRUE.
#' @param validation (Only for "plsda" and "splsda" methods) Validation method. Options are "Mfold" and "loo".
#' @param folds (Only for "plsda" and "splsda" methods) Numeric. Number of folds for Mfold validation method (default is 5). If the validation method is loo, this value will become to 1.
#' @param nrepeat (Only for "plsda" and "splsda" methods) Numeric. Number of iterations for the validation method selected.
#' @param vip (Only for "plsda" method) Numeric indicating VIP cutoff to select features that will be displayed in vip plot.
#' @param num_features (Only for "splsda" method) Numeric. Number of variables selected to discriminate groups.
#' @param legend_position Character indicating the legend position. Options are "none", "top", "bottom", "left", and "right".
#'
#' @export
#'
#' @return A list with all results for multivariate statistical analysis including plots and data frames.
#' @author Pol Castellano-Escuder
#'
#' @import ggplot2
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr select arrange desc filter
#' @importFrom tidyr pivot_longer
#' @importFrom magrittr %>%
#' @importFrom mixOmics pca plsda perf vip tune.splsda splsda selectVar
#' @importFrom MSnbase pData exprs
#' 
#' @examples 
#' data("st000336")
#' 
#' # PCA
#' st000336 %>% 
#'   PomaImpute() %>%
#'   PomaNorm() %>%
#'   PomaOutliers() %>%
#'   PomaMultivariate(method = "pca")
#' 
#' # PLSDA
#' st000336 %>% 
#'   PomaImpute() %>%
#'   PomaNorm() %>%
#'   PomaOutliers() %>%
#'   PomaMultivariate(method = "plsda", vip = 1)
#' 
#' # sPLSDA
#' st000336 %>% 
#'   PomaImpute() %>%
#'   PomaNorm() %>%
#'   PomaOutliers() %>%
#'   PomaMultivariate(method = "splsda")
PomaMultivariate <- function(data,
                             method = "pca",
                             components = 5,
                             center = FALSE,
                             scale = FALSE,
                             labels = FALSE,
                             load_length = 1,
                             ellipse = TRUE,
                             validation = "Mfold",
                             folds = 5,
                             nrepeat = 10,
                             vip = 1.5,
                             num_features = 10,
                             legend_position = "bottom"){

  if (missing(data)) {
    stop("data argument is empty!")
  }
  if(!is(data[1], "SummarizedExperiment")){
    stop("data is not a SummarizedExperiment object. \nSee POMA::PomaSummarizedExperiment or SummarizedExperiment::SummarizedExperiment")
  }
  if (missing(method)) {
    stop("Select a method!")
  }
  if (!(method %in% c("pca", "plsda", "splsda"))) {
    stop("Incorrect value for method argument!")
  }
  if (missing(validation)) {
    if (method %in% c("plsda", "splsda")){
      message("validation argument is empty! Mfold will be used")
    }
  }
  if (!(validation %in% c("Mfold", "loo"))) {
    stop("Incorrect validation method! Please choose 'Mfold' or 'loo'")
  }
  if (load_length > 2 | load_length < 1) {
    stop("load_length should be a number between 1 and 2...")
  }
  if(!(legend_position %in% c("none", "top", "bottom", "left", "right"))) {
    stop("Incorrect value for legend_position argument!")
  }

  df <- t(MSnbase::exprs(data))

  if(method == "pca"){

    X <- as.matrix(df)
    Y <- as.factor(MSnbase::pData(data)[,1])
    pca_res <- mixOmics::pca(X, ncomp = components, center = center, scale = scale)

    PCi <- data.frame(pca_res$variates$X, Groups = Y) %>% 
      rownames_to_column("ID")

    scoresplot <- ggplot(PCi, aes(x = PC1, y = PC2, color = Groups, shape = Groups, label = ID)) +
      {if(!labels)geom_point(size = 3, alpha = 0.8)} +
      xlab(paste0("PC1 (", round(100*(pca_res$prop_expl_var$X)[1], 2), "%)")) +
      ylab(paste0("PC2 (", round(100*(pca_res$prop_expl_var$X)[2], 2), "%)")) +
      {if(ellipse)stat_ellipse(type = "norm")} +
      {if(labels)geom_text(aes(label = ID), show.legend = TRUE)} +
      theme_bw() +
      theme(legend.title = element_blank(),
            legend.position = legend_position) +
      scale_colour_viridis_d()

    ##

    eigenvalues <- data.frame(Percent_Variance_Explained = round(pca_res$prop_expl_var$X*100, 3)) %>% 
      rownames_to_column("Principal_Component")

    screeplot <- ggplot(eigenvalues, aes(x = reorder(Principal_Component, -Percent_Variance_Explained), y = Percent_Variance_Explained, fill = Percent_Variance_Explained)) +
      geom_bar(stat = "identity") +
      xlab("Principal Component") +
      ylab("% Variance Explained") +
      theme_bw() +
      theme(legend.title = element_blank(),
            legend.position = legend_position) +
      scale_fill_continuous(type = "viridis")

    ####

    score_data <- PCi %>% dplyr::select(-Groups, -ID) %>% round(4)

    ####

    pca_res2 <- mixOmics::pca(X, ncomp = components, center = TRUE, scale = TRUE)

    PCi2 <- data.frame(pca_res2$variates$X, Groups = Y)

    lam <- (pca_res2$sdev[1:2] * sqrt(nrow(PCi2)))^load_length
    len <- t(t(pca_res2$loadings$X[, 1:2]) * lam)*0.8
    PCAloadings <- data.frame(pca_res2$loadings$X, to_x = len[,1], to_y = len[,2])

    biplot <- ggplot(PCi2, aes(x = PC1, y = PC2, color = Groups))+
      geom_point(size = 3, alpha = 0.8) +
      xlab(paste0("PC1 (", round(100*(pca_res2$prop_expl_var$X)[1], 2), "%)")) +
      ylab(paste0("PC2 (", round(100*(pca_res2$prop_expl_var$X)[2], 2), "%)")) +
      theme_bw() +
      {if(ellipse)stat_ellipse(type = "norm")} +
      geom_segment(data = PCAloadings,
                   aes(x = 0, y = 0,
                       xend = to_x,
                       yend = to_y),
                   arrow = arrow(length = unit(1/2, "picas")), color = "grey19") +
      annotate("text", x = PCAloadings$to_x,
               y = PCAloadings$to_y,
               label = rownames(PCAloadings), size = 4) +
      theme(legend.title = element_blank(),
            legend.position = legend_position) +
      scale_colour_viridis_d()

    return(list(screeplot = screeplot, scoresplot = scoresplot,
                score_data = score_data, eigenvalues = eigenvalues, biplot = biplot))

  }

  else if (method == "plsda"){

    X <- as.matrix(df)
    Y <- as.factor(MSnbase::pData(data)[,1])

    plsda_res <- mixOmics::plsda(X, Y, ncomp = components)

    PLSDAi <- data.frame(plsda_res$variates$X, Groups = Y) %>% 
      rownames_to_column("ID")

    scoresplot <- ggplot(PLSDAi, aes(x = comp1, y = comp2, color = Groups, shape = Groups, label = ID))+
      {if(!labels)geom_point(size = 3, alpha = 0.8)} +
      xlab("Component 1") +
      ylab("Component 2") +
      {if(ellipse)stat_ellipse(type = "norm")} +
      {if(labels)geom_text(aes(label = ID), show.legend = TRUE)} +
      theme_bw() +
      theme(legend.title = element_blank(),
            legend.position = legend_position) +
      scale_colour_viridis_d()

    #####

    perf_plsda <- mixOmics::perf(plsda_res, validation = validation, folds = folds,
                                 progressBar = TRUE, auc = TRUE, nrepeat = nrepeat)

    overall <- data.frame(perf_plsda$error.rate[1]) %>%
      round(4) %>%
      rownames_to_column("Component") %>%
      pivot_longer(cols = -Component)

    ber <- data.frame(perf_plsda$error.rate[2]) %>%
      round(4) %>%
      rownames_to_column("Component") %>%
      pivot_longer(cols = -Component)

    errors_plsda <- rbind(ber, overall)

    errors_plsda_plot <- ggplot(data = errors_plsda, aes(x = Component, y = value, group = name)) +
      geom_line(aes(color = name)) +
      geom_point(aes(color = name)) +
      ylab("Error") +
      theme_bw() +
      theme(legend.title = element_blank(),
            legend.position = legend_position) +
      scale_colour_viridis_d()

    ####

    plsda_vip <- data.frame(mixOmics::vip(plsda_res)) %>%
      rownames_to_column("feature") %>%
      arrange(desc(comp1))

    plsda_vip_top <- plsda_vip %>%
      filter(comp1 >= vip) %>%
      mutate(feature = factor(feature, levels = feature[order(comp1)]))

    vip_plsda_plot <- ggplot(plsda_vip_top, aes(x = reorder(feature, comp1), y = comp1, fill = comp1)) +
      geom_bar(stat = "identity") +
      coord_flip() +
      ylab("VIP") +
      xlab("") +
      # {if(nrow(plsda_vip_top) <= 10)geom_label(data = plsda_vip_top, aes(label = round(comp1, 2)))} +
      theme_bw() +
      theme(legend.title = element_blank(),
            legend.position = legend_position) +
      scale_fill_continuous(type = "viridis")

    ####

    scores_plsda <- PLSDAi %>% dplyr::select(-Groups, -ID) %>% round(4)

    return(list(scoresplot = scoresplot, errors_plsda = errors_plsda,
                errors_plsda_plot = errors_plsda_plot, plsda_vip_table = plsda_vip,
                vip_plsda_plot = vip_plsda_plot, score_data = scores_plsda))
  }

  else if (method == "splsda"){

    X <- as.matrix(df)
    Y <- as.factor(MSnbase::pData(data)[,1])

    list_keepX <- c(1:num_features)

    tune_splsda <- mixOmics::tune.splsda(X, Y, ncomp = components, validation = validation, folds = folds,
                                         progressBar = TRUE, dist = 'max.dist', measure = "BER",
                                         test.keepX = list_keepX, nrepeat = nrepeat) # cpus = 4

    error <- tune_splsda$error.rate
    ncomp <- tune_splsda$choice.ncomp$ncomp # optimal number of components based on t-tests
    select_keepX <- tune_splsda$choice.keepX[1:ncomp]  # optimal number of variables to select

    errors_splsda <- data.frame(tune_splsda$error.rate) %>% 
      round(4) %>%
      rownames_to_column("features") %>%
      pivot_longer(cols = -features)

    errors_sd <- data.frame(tune_splsda$error.rate.sd) %>% 
      rownames_to_column("features_sd") %>%
      pivot_longer(cols = -features_sd)

    errors_splsda <- cbind(errors_splsda, sd = errors_sd$value)

    bal_error_rate <- ggplot(data = errors_splsda, aes(x = features, y = value, group = name)) +
      geom_line(aes(color = name)) +
      geom_point(aes(color = name)) +
      geom_errorbar(aes(ymin = value-sd, ymax = value+sd, color = name), width=.1) +
      theme_bw() +
      xlab("Number of features") +
      ylab("Error") +
      theme(legend.title = element_blank(),
            legend.position = legend_position) +
      scale_colour_viridis_d()

    ####

    if (ncomp == 1){
      ncompX <- 2
    }else{
      ncompX <- ncomp}

    ####

    res_splsda <- mixOmics::splsda(X, Y, ncomp = ncompX, keepX = select_keepX)

    SPLSDAi <- data.frame(res_splsda$variates$X, Groups = Y) %>% rownames_to_column("ID")

    splsda_scores_plot <- ggplot(SPLSDAi, aes(x = comp1, y = comp2, color = Groups, shape = Groups, label = ID)) +
      {if(!labels)geom_point(size = 3, alpha = 0.8)} +
      xlab("Component 1") +
      ylab("Component 2") +
      {if(ellipse)stat_ellipse(type = "norm")} +
      {if(labels)geom_text(aes(label = ID), show.legend = TRUE)} +
      theme_bw() +
      theme(legend.title = element_blank(),
            legend.position = legend_position) +
      scale_colour_viridis_d()

    scores_splsda <- SPLSDAi %>% dplyr::select(-Groups, -ID) %>% round(4)

    selected_variables <- mixOmics::selectVar(res_splsda, comp = 1)
    selected_variables <- round(selected_variables$value, 4)
    selected_variables <- data.frame(Feature = rownames(selected_variables), Value = selected_variables$value.var)

    ####

    return(list(ncomp = ncomp, select_keepX = select_keepX, errors_splsda = errors_splsda,
                scoresplot = splsda_scores_plot, bal_error_rate = bal_error_rate,
                score_data = scores_splsda, selected_variables = selected_variables))

  }

}

