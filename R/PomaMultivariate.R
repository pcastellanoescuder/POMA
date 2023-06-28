
#' Multivariate Statistical Methods for Mass Spectrometry Data
#'
#' @description PomaMultivariate() allows users to perform different multivariate statistical analysis on MS data.
#'
#' @param data A SummarizedExperiment object.
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
#' @return A list with all results for multivariate statistical analysis including plots and tables.
#' @author Pol Castellano-Escuder
#'
#' @importFrom magrittr %>%
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
  if(!is(data, "SummarizedExperiment")){
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

  df <- t(SummarizedExperiment::assay(data))

  if (method == "plsda"){

    X <- as.matrix(df)
    Y <- as.factor(SummarizedExperiment::colData(data)[,1])

    plsda_res <- mixOmics::plsda(X, Y, ncomp = components)

    PLSDAi <- data.frame(plsda_res$variates$X, Groups = Y) %>% 
      tibble::rownames_to_column("ID")

    scoresplot <- ggplot2::ggplot(PLSDAi, ggplot2::aes(x = comp1, y = comp2, color = Groups, shape = Groups, label = ID))+
      {if(!labels)ggplot2::geom_point(size = 2, alpha = 0.9)} +
      ggplot2::labs(x = "Component 1",
                    y = "Component 2") +
      {if(ellipse)ggplot2::stat_ellipse(type = "norm")} +
      {if(labels)ggplot2::geom_text(ggplot2::aes(label = ID), show.legend = TRUE)} +
      ggplot2::theme_bw() +
      ggplot2::theme(legend.title = ggplot2::element_blank(),
                     legend.position = legend_position) +
      ggplot2::scale_color_viridis_d(option = "plasma", end = 0.8)

    perf_plsda <- mixOmics::perf(plsda_res, validation = validation, folds = folds,
                                 progressBar = TRUE, auc = TRUE, nrepeat = nrepeat)

    overall <- data.frame(perf_plsda$error.rate[1]) %>%
      round(4) %>%
      tibble::rownames_to_column("component") %>%
      tidyr::pivot_longer(cols = -component) %>% 
      dplyr::as_tibble()

    ber <- data.frame(perf_plsda$error.rate[2]) %>%
      round(4) %>%
      tibble::rownames_to_column("component") %>%
      tidyr::pivot_longer(cols = -component) %>% 
      dplyr::as_tibble()

    errors_plsda <- rbind(ber, overall)

    errors_plsda_plot <- ggplot2::ggplot(data = errors_plsda, ggplot2::aes(x = component, y = value, group = name)) +
      ggplot2::geom_line(ggplot2::aes(color = name)) +
      ggplot2::geom_point(ggplot2::aes(color = name)) +
      ggplot2::labs(x = "Component",
                    y = "Error") +
      ggplot2::theme_bw() +
      ggplot2::theme(legend.title = ggplot2::element_blank(),
                     legend.position = legend_position) +
      ggplot2::scale_fill_viridis_d(option = "plasma", end = 0.8)

    plsda_vip <- data.frame(mixOmics::vip(plsda_res)) %>%
      tibble::rownames_to_column("feature") %>%
      dplyr::arrange(dplyr::desc(comp1)) %>% 
      dplyr::as_tibble()

    plsda_vip_top <- plsda_vip %>%
      dplyr::filter(comp1 >= vip) %>%
      dplyr::mutate(feature = factor(feature, levels = feature[order(comp1)]))

    vip_plsda_plot <- ggplot2::ggplot(plsda_vip_top, ggplot2::aes(x = reorder(feature, comp1), y = comp1, fill = comp1)) +
      ggplot2::geom_bar(stat = "identity") +
      ggplot2::coord_flip() +
      ggplot2::labs(x = NULL,
                    y = "VIP") +
      ggplot2::theme_bw() +
      ggplot2::theme(legend.title = ggplot2::element_blank(),
                     legend.position = legend_position) +
      ggplot2::scale_fill_viridis_c(option = "plasma", end = 0.8)
    
    scores_plsda <- PLSDAi %>% 
      dplyr::select(-Groups, -ID) %>%
      dplyr::as_tibble()

    return(list(scoresplot = scoresplot, 
                errors_plsda = errors_plsda,
                errors_plsda_plot = errors_plsda_plot,
                vip_plsda = plsda_vip,
                vip_plsda_plot = vip_plsda_plot, 
                scores = scores_plsda))
  }

  else if (method == "splsda"){

    X <- as.matrix(df)
    Y <- as.factor(SummarizedExperiment::colData(data)[,1])

    list_keepX <- c(1:num_features)

    tune_splsda <- mixOmics::tune.splsda(X, Y, ncomp = components, validation = validation, folds = folds,
                                         progressBar = TRUE, dist = 'max.dist', measure = "BER",
                                         test.keepX = list_keepX, nrepeat = nrepeat) # cpus = 4

    error <- tune_splsda$error.rate
    ncomp <- tune_splsda$choice.ncomp$ncomp # optimal number of components based on t-tests
    select_keepX <- tune_splsda$choice.keepX[1:ncomp]  # optimal number of variables to select

    errors_splsda <- data.frame(tune_splsda$error.rate) %>% 
      tibble::rownames_to_column("feature") %>%
      tidyr::pivot_longer(cols = -feature)

    errors_sd <- data.frame(tune_splsda$error.rate.sd) %>% 
      tibble::rownames_to_column("feature_sd") %>%
      tidyr::pivot_longer(cols = -feature_sd)

    errors_splsda <- cbind(errors_splsda, sd = errors_sd$value) %>% 
      as.data.frame() %>% 
      dplyr::as_tibble()

    bal_error_rate <- ggplot2::ggplot(data = errors_splsda, ggplot2::aes(x = feature, y = value, group = name)) +
      ggplot2::geom_line(ggplot2::aes(color = name)) +
      ggplot2::geom_point(ggplot2::aes(color = name)) +
      ggplot2::geom_errorbar(ggplot2::aes(ymin = value - sd, ymax = value + sd, color = name), width = 0.1) +
      ggplot2::theme_bw() +
      ggplot2::labs(x = "Number of features",
                    y = "Error") +
      ggplot2::theme(legend.title = ggplot2::element_blank(),
                     legend.position = legend_position) +
      ggplot2::scale_color_viridis_d(option = "plasma", end = 0.8)

    if (ncomp == 1) {
      ncompX <- 2
      } else {
        ncompX <- ncomp
      }

    res_splsda <- mixOmics::splsda(X, Y, ncomp = ncompX, keepX = select_keepX)

    SPLSDAi <- data.frame(res_splsda$variates$X, Groups = Y) %>% 
      tibble::rownames_to_column("ID")

    splsda_scores_plot <- ggplot2::ggplot(SPLSDAi, ggplot2::aes(x = comp1, y = comp2, color = Groups, shape = Groups, label = ID)) +
      {if(!labels)ggplot2::geom_point(size = 2, alpha = 0.9)} +
      ggplot2::labs(x = "Component 1",
                    y = "Component 2") +
      {if(ellipse)ggplot2::stat_ellipse(type = "norm")} +
      {if(labels)ggplot2::geom_text(ggplot2::aes(label = ID), show.legend = TRUE)} +
      ggplot2::theme_bw() +
      ggplot2::theme(legend.title = ggplot2::element_blank(),
                     legend.position = legend_position) +
      ggplot2::scale_color_viridis_d(option = "plasma", end = 0.8)

    scores_splsda <- SPLSDAi %>% 
      dplyr::select(-Groups, -ID) %>%
      dplyr::as_tibble()

    selected_variables <- mixOmics::selectVar(res_splsda, comp = 1)$value %>% 
      dplyr::mutate(value.var  = round(value.var, 4)) %>% 
      dplyr::rename(value = value.var) %>% 
      tibble::rownames_to_column("feature") %>% 
      dplyr::as_tibble()

    return(list(ncomp = ncomp, 
                select_keepX = select_keepX, 
                scoresplot = splsda_scores_plot, 
                errors_splsda = errors_splsda,
                errors_splsda_plot = bal_error_rate,
                scores = scores_splsda, 
                selected_variables = selected_variables))

  }

}

