
#' Partial Least Squares Methods
#'
#' @description `PomaPLS` performs Partial Least Squares (PLS) regression, Partial Least Squares Discriminant Analysis (PLS-DA) to classify samples, and Sparse Partial Least Squares Discriminant Analysis (sPLS-DA) to classify samples (supervised analysis) and select variables.
#'
#' @param data A `SummarizedExperiment` object.
#' @param method Character. PLS method. Options include "pls", "plsda", and "splsda".
#' @param y Character. Indicates the name of `colData` columns to be used as dependent variable. If it's set to NULL, the first variable in `colData` will be used as the dependent variable.
#' @param ncomp Numeric. Number of components in the model. Default is 5.
#' @param labels Logical. Indicates if sample names should be displayed.
#' @param ellipse Logical. Indicates whether a 95 percent confidence interval ellipse should be displayed. Default is TRUE.
#' @param cross_validation Logical. Indicates if cross-validation should be performed for PLS-DA ("plsda") and sPLS-DA ("splsda") methods. Default is FALSE.
#' @param validation Character. (Only for "plsda" and "splsda" methods). Indicates the cross-validation method. Options are "Mfold" and "loo" (Leave-One-Out).
#' @param folds Numeric. (Only for "plsda" and "splsda" methods). Number of folds for "Mfold" cross-validation method (default is 5). If the validation method is "loo", this value is set to 1.
#' @param nrepeat Numeric. (Only for "plsda" and "splsda" methods). Number of times the cross-validation process is repeated.
#' @param vip Numeric. (Only for "plsda" method). Indicates the variable importance in the projection (VIP) cutoff.
#' @param num_features Numeric. (Only for "splsda" method). Number of features to discriminate groups.
#' @param theme_params List. Indicates `theme_poma` parameters.
#' 
#' @export
#'
#' @return A `list` with results including plots and tables.
#' @author Pol Castellano-Escuder
#'
#' @importFrom magrittr %>%
#' 
#' @examples 
#' data("st000284")
#' 
#' # PLS
#' st000284 %>% 
#'   PomaNorm() %>%
#'   PomaPLS(method = "pls")
#'   
#' data("st000336")
#' 
#' # PLSDA
#' st000336 %>% 
#'   PomaImpute() %>%
#'   PomaNorm() %>%
#'   PomaPLS(method = "plsda")
#'   
#' # PLSDA with Cross-Validation
#' st000336 %>% 
#'   PomaImpute() %>%
#'   PomaNorm() %>%
#'   PomaPLS(method = "plsda", cross_validation = TRUE)
#'   
#' # sPLSDA
#' st000336 %>% 
#'   PomaImpute() %>%
#'   PomaNorm() %>%
#'   PomaPLS(method = "splsda")
#'   
#' # sPLSDA with Cross-Validation
#' st000336 %>% 
#'   PomaImpute() %>%
#'   PomaNorm() %>%
#'   PomaPLS(method = "splsda", ncomp = 3, cross_validation = TRUE)
PomaPLS <- function(data,
                    method = "pls",
                    y = NULL,
                    ncomp = 5,
                    labels = FALSE,
                    ellipse = TRUE,
                    cross_validation = FALSE,
                    validation = "Mfold",
                    folds = 5,
                    nrepeat = 10,
                    vip = 1,
                    num_features = 10,
                    theme_params = list(),
                    ...) {

  if(!is(data, "SummarizedExperiment")){
    stop("data is not a SummarizedExperiment object. \nSee POMA::PomaCreateObject or SummarizedExperiment::SummarizedExperiment")
  }
  if (!(method %in% c("pls", "plsda", "splsda"))) {
    stop("Incorrect value for method argument")
  }
  if (!(validation %in% c("Mfold", "loo"))) {
    stop("Incorrect value for validation argument. Options are: 'Mfold' and 'loo'")
  }

  to_pls <- t(SummarizedExperiment::assay(data))

  if (method == "pls") {
    dependent_variable <- SummarizedExperiment::colData(data) %>% 
      as.data.frame() %>% 
      dplyr::select_if(is.numeric)
    
    if (ncol(dependent_variable) == 0) {
      stop("No numeric variables to be used as dependent variable in metadata file")
    }
    if (is.null(y)) {
      y <- colnames(dependent_variable)[1]
    }
    
    dependent_variable <- dependent_variable %>% 
      dplyr::select(dplyr::all_of(y)) %>% 
      as.matrix()
    
    pls_res <- mixOmics::pls(X = to_pls,
                             Y = dependent_variable,
                             ncomp = ncomp)
    
    # factors
    pls_res_df <- data.frame(y = as.numeric(dependent_variable), pls_res$variates$X) %>% 
      tibble::rownames_to_column("sample_id") %>% 
      dplyr::as_tibble()
    
    # factors plot
    factors_plot <- ggplot2::ggplot(pls_res_df, ggplot2::aes(x = comp1, y = comp2))+
      {if(!labels)ggplot2::geom_point(ggplot2::aes(fill = y), size = 3, alpha = 0.6, pch = 21)} +
      ggplot2::labs(x = "Component 1",
                    y = "Component 2",
                    fill = "Dependent\nvariable (Y)",
                    color = "Dependent\nvariable (Y)") +
      {if(ellipse)ggplot2::stat_ellipse(type = "norm", show.legend = FALSE)} +
      {if(labels)ggplot2::geom_text(ggplot2::aes(color = y, label = sample_id))} +
      do.call(theme_poma, theme_params) +
      scale_fill_poma_c() +
      scale_color_poma_c()
    
    # loadings
    pls_loadings_df <- data.frame(pls_res$loadings$X) %>% 
      tibble::rownames_to_column("feature") %>% 
      dplyr::as_tibble()
    
    # loadings plot
    loadings_plot <- pls_loadings_df %>%
      dplyr::arrange(dplyr::desc(abs(comp1))) %>%
      dplyr::select(feature, comp1:paste0("comp", ncomp)) %>% 
      dplyr::slice(1L:10L) %>%
      tidyr::pivot_longer(cols = -feature) %>% 
      ggplot2::ggplot(ggplot2::aes(x = reorder(feature, value), 
                                   y = value,
                                   fill = name)) +
      ggplot2::geom_col(position = "dodge2") +
      ggplot2::labs(x = NULL,
                    y = "Loadings",
                    fill = NULL) +
      theme_poma(axis_x_rotate = TRUE) +
      scale_fill_poma_d()
    
    return(list(factors = pls_res_df,
                factors_plot = factors_plot,
                loadings = pls_loadings_df,
                loadings_plot = loadings_plot)
           )
  }
  
  else if (method == "plsda") {
    dependent_variable <- SummarizedExperiment::colData(data) %>% 
      as.data.frame() %>% 
      dplyr::select_if(is.factor)
    
    if (ncol(dependent_variable) == 0) {
      stop("No factor variables to be used as dependent variable in metadata file")
    }
    if (is.null(y)) {
      y <- colnames(dependent_variable)[1]
    }
    
    dependent_variable <- dependent_variable %>% 
      dplyr::select(dplyr::all_of(y[1])) %>% 
      dplyr::pull(1)

    plsda_res <- mixOmics::plsda(X = to_pls,
                                 Y = dependent_variable,
                                 ncomp = ncomp)
    
    # factors
    plsda_res_df <- data.frame(y = dependent_variable, plsda_res$variates$X) %>% 
      tibble::rownames_to_column("sample_id") %>% 
      dplyr::as_tibble()
   
    # factors plot
    factors_plot <- ggplot2::ggplot(plsda_res_df, ggplot2::aes(x = comp1, y = comp2))+
      {if(!labels)ggplot2::geom_point(ggplot2::aes(fill = y), size = 3, alpha = 0.6, pch = 21)} +
      ggplot2::labs(x = "Component 1",
                    y = "Component 2",
                    fill = "Dependent\nvariable (Y)",
                    color = "Dependent\nvariable (Y)") +
      {if(ellipse)ggplot2::stat_ellipse(ggplot2::aes(color = y), type = "norm", show.legend = FALSE)} +
      {if(labels)ggplot2::geom_text(ggplot2::aes(color = y, label = sample_id))} +
      do.call(theme_poma, theme_params) +
      scale_fill_poma_d() +
      scale_color_poma_d()
    
    # VIP
    plsda_vip_df <- data.frame(mixOmics::vip(plsda_res)) %>%
      tibble::rownames_to_column("feature") %>%
      dplyr::arrange(dplyr::desc(comp1)) %>% 
      dplyr::as_tibble()

    # VIP plot
    plsda_vip_top <- plsda_vip_df %>%
      dplyr::filter(comp1 >= vip)

    vip_plot <- ggplot2::ggplot(plsda_vip_top, ggplot2::aes(x = comp1, y = reorder(feature, comp1), fill = comp1)) +
      ggplot2::geom_col() +
      ggplot2::labs(x = "Variable Importance\nin the Projection (VIP)",
                    y = NULL,
                    fill = "Component 1\nVIP") +
      theme_poma(legend_position = "right") +
      scale_fill_poma_c()

    # cross-validation
    if (cross_validation) {
      perf_plsda <- mixOmics::perf(plsda_res, validation = validation, folds = folds,
                                   progressBar = TRUE, auc = TRUE, nrepeat = nrepeat)

      ber_errors <- data.frame(perf_plsda$error.rate$BER) %>%
        janitor::clean_names() %>% 
        tibble::rownames_to_column("component") %>%
        tidyr::pivot_longer(cols = -component) %>%
        dplyr::mutate(error = "BER") %>% 
        dplyr::as_tibble()

      errors_plsda <- data.frame(perf_plsda$error.rate$overall) %>%
        janitor::clean_names() %>% 
        tibble::rownames_to_column("component") %>%
        tidyr::pivot_longer(cols = -component) %>%
        dplyr::mutate(error = "Overall") %>% 
        dplyr::bind_rows(ber_errors)
      
      # errors plot
      errors_plsda_plot <- ggplot2::ggplot(data = errors_plsda, ggplot2::aes(x = component, y = value)) +
        ggplot2::geom_line(ggplot2::aes(color = name, group = name), show.legend = FALSE) +
        ggplot2::geom_point(ggplot2::aes(fill = name), size = 3, alpha = 0.6, pch = 21) +
        ggplot2::labs(x = "Component",
                      y = "Error",
                      fill = "Error\nType") +
        ggplot2::facet_wrap(~error) +
        do.call(theme_poma, theme_params) +
        scale_color_poma_d() +
        scale_fill_poma_d()
      
      # errors
      errors <- errors_plsda %>% 
        dplyr::mutate(component = paste0(component, "_", error)) %>% 
        dplyr::select(-error) %>% 
        tidyr::pivot_wider(id_cols = component, names_from = name, values_from = value) %>% 
        dplyr::as_tibble()
      
      return(list(factors = plsda_res_df,
                  factors_plot = factors_plot,
                  vip_values = plsda_vip_df,
                  vip_plot = vip_plot,
                  errors = errors_plsda,
                  errors_plot = errors_plsda_plot)
             )
      
    } else {
      return(list(factors = plsda_res_df,
                  factors_plot = factors_plot,
                  vip_values = plsda_vip_df,
                  vip_plot = vip_plot)
             ) 
    }
  }

  else if (method == "splsda") {
    dependent_variable <- SummarizedExperiment::colData(data) %>% 
      as.data.frame() %>% 
      dplyr::select_if(is.factor)
    
    if (ncol(dependent_variable) == 0) {
      stop("No factor variables to be used as dependent variable in metadata file")
    }
    if (is.null(y)) {
      y <- colnames(dependent_variable)[1]
    }
    
    dependent_variable <- dependent_variable %>% 
      dplyr::select(dplyr::all_of(y[1])) %>% 
      dplyr::pull(1)
    
    splsda_res <- mixOmics::splsda(X = to_pls,
                                   Y = dependent_variable,
                                   ncomp = ncomp,
                                   keepX = num_features)

    # factors
    splsda_res_df <- data.frame(y = dependent_variable, splsda_res$variates$X) %>% 
      tibble::rownames_to_column("sample_id") %>% 
      dplyr::as_tibble()
    
    # factors plot
    factors_plot <- ggplot2::ggplot(splsda_res_df, ggplot2::aes(x = comp1, y = comp2))+
      {if(!labels)ggplot2::geom_point(ggplot2::aes(fill = y), size = 3, alpha = 0.6, pch = 21)} +
      ggplot2::labs(x = "Component 1",
                    y = "Component 2",
                    fill = "Dependent\nvariable (Y)",
                    color = "Dependent\nvariable (Y)") +
      {if(ellipse)ggplot2::stat_ellipse(ggplot2::aes(color = y), type = "norm", show.legend = FALSE)} +
      {if(labels)ggplot2::geom_text(ggplot2::aes(color = y, label = sample_id))} +
      do.call(theme_poma, theme_params) +
      scale_fill_poma_d() +
      scale_color_poma_d()
    
    # selected features
    selected_features <- mixOmics::selectVar(splsda_res, comp = 1)$value %>% 
      dplyr::rename(value = 1) %>% 
      tibble::rownames_to_column("feature") %>% 
      dplyr::as_tibble()
    
    # selected features plot
    selected_features_plot <- ggplot2::ggplot(selected_features, ggplot2::aes(x = value, y = reorder(feature, value), fill = value)) +
      ggplot2::geom_col() +
      ggplot2::labs(x = "Value",
                    y = NULL,
                    fill = "Component 1\nValue") +
      theme_poma(legend_position = "right") +
      scale_fill_poma_c()
    
    # cross-validation
    if (cross_validation) {
      tune_splsda <- mixOmics::tune.splsda(X = to_pls,
                                           Y = dependent_variable, 
                                           ncomp = ncomp, validation = validation, folds = folds,
                                           progressBar = TRUE, dist = 'max.dist', measure = "BER",
                                           test.keepX = c(1:num_features), nrepeat = nrepeat)

      opt_ncomp <- tune_splsda$choice.ncomp$ncomp # optimal number of components based on t-tests
      select_keepX <- tune_splsda$choice.keepX[1:ncomp]  # optimal number of variables to select

      errors_splsda_sd <- data.frame(tune_splsda$error.rate.sd) %>%
        tibble::rownames_to_column("feature_sd") %>%
        tidyr::pivot_longer(cols = -feature_sd) %>% 
        dplyr::select(sd = value)
      
      errors_splsda <- data.frame(tune_splsda$error.rate) %>%
        tibble::rownames_to_column("feature") %>%
        tidyr::pivot_longer(cols = -feature) %>% 
        dplyr::bind_cols(errors_splsda_sd) %>% 
        dplyr::as_tibble()
      
      # errors plot
      errors_splsda_plot <- ggplot2::ggplot(errors_splsda, ggplot2::aes(x = as.factor(as.numeric(feature)), y = value)) +
        ggplot2::geom_line(ggplot2::aes(color = name, group = name), show.legend = FALSE) +
        ggplot2::geom_point(ggplot2::aes(fill = name), size = 3, alpha = 0.6, pch = 21) +
        ggplot2::geom_errorbar(ggplot2::aes(ymin = value - sd, ymax = value + sd, color = name), width = 0.1, show.legend = FALSE) +
        ggplot2::labs(x = "# Features",
                      y = "Error",
                      fill = "Component") +
        do.call(theme_poma, theme_params) +
        scale_color_poma_d() +
        scale_fill_poma_d()
      
      # errors
      errors <- errors_splsda %>% 
        dplyr::select(feature, value, name) %>% 
        tidyr::pivot_wider(id_cols = feature, names_from = name, values_from = value) %>% 
        dplyr::as_tibble()
      
      return(list(factors = splsda_res_df,
                  factors_plot = factors_plot,
                  selected_features = selected_features,
                  selected_features_plot = selected_features_plot,
                  errors = errors_splsda,
                  errors_plot = errors_splsda_plot,
                  optimal_components = opt_ncomp,
                  optimal_features = select_keepX)
             )
      
    } else {
      return(list(factors = splsda_res_df,
                  factors_plot = factors_plot,
                  selected_features = selected_features,
                  selected_features_plot = selected_features_plot)
             ) 
    }
  }
}

