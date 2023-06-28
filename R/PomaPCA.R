
#' Principal Components Analysis
#'
#' @description `PomaPCA` performs a principal components analysis on the given `SummarizedExperiment` object.
#'
#' @param data A `SummarizedExperiment` object.
#' @param center Logical. Indicates whether the variables should be shifted to be zero centered. Default is TRUE.
#' @param scale Logical. Indicates whether the variables should be scaled to have unit variance before the analysis takes place. Default is TRUE.
#' @param labels Logical. Indicates if sample names should be displayed.
#' @param ellipse Logical. Indicates whether a 95 percent confidence interval ellipse should be displayed in score plot and biplot. Default is FALSE.
#' @param load_length Numeric. Indicates the length of biplot loading arrows. Value between 1 and 2. Default is 1.
#' 
#' @export
#'
#' @return A `list` with results including plots and tables.
#' @author Pol Castellano-Escuder
#'
#' @importFrom magrittr %>%
#' 
#' @examples 
#' data("st000336")
#' 
#' st000336 %>% 
#'   PomaImpute() %>%
#'   PomaNorm() %>%
#'   PomaPCA()
PomaPCA <- function(data,
                    center = TRUE,
                    scale = TRUE,
                    labels = FALSE,
                    ellipse = FALSE,
                    load_length = 1,
                    ...) {
  
  if(!is(data, "SummarizedExperiment")){
    stop("data is not a SummarizedExperiment object. \nSee POMA::PomaCreateObject or SummarizedExperiment::SummarizedExperiment")
  }
  
  to_pca <- t(SummarizedExperiment::assay(data))
  pca_res <- prcomp(to_pca, center = center, scale. = scale)
  
  if (ncol(SummarizedExperiment::colData(data)) != 0) {
    group_factor <- SummarizedExperiment::colData(data)[,1]
    pca_res_df <- data.frame(sample_id = rownames(SummarizedExperiment::colData(data)),
                             group = group_factor, 
                             pca_res$x)
  } else {
    pca_res_df <- data.frame(sample_id = colnames(SummarizedExperiment::assay(data)),
                             pca_res$x)
  }
  
  # scores plot
  scores_plot <- ggplot2::ggplot(pca_res_df, ggplot2::aes(x = PC1, y = PC2)) +
    {if(ncol(SummarizedExperiment::colData(data)) != 0)ggplot2::geom_point(ggplot2::aes(fill = group), pch = 21, size = 3, alpha = 0.6)} +
    {if(ncol(SummarizedExperiment::colData(data)) == 0)ggplot2::geom_point(pch = 21, size = 3, alpha = 0.6)} +
    ggplot2::labs(x = paste0("PC1 (", round(100*(((pca_res$sdev^2)[1]) / sum(pca_res$sdev^2)), 2), "%)"),
                  y = paste0("PC2 (", round(100*(((pca_res$sdev^2)[2]) / sum(pca_res$sdev^2)), 2), "%)"),
                  fill = NULL,
                  color = NULL) +
    {if(ellipse)ggplot2::stat_ellipse(ggplot2::aes(color = group), type = "norm")} +
    {if(labels)ggplot2::geom_text(ggplot2::aes(label = sample_id), show.legend = FALSE)} +
    theme_poma() +
    scale_fill_poma_d() +
    scale_color_poma_d()
  
  # scree plot
  eigenvalues <- data.frame(comp = paste0("PC", 1:ncol(pca_res$x)),
                            var_exp = round(100*(((pca_res$sdev^2)) / sum(pca_res$sdev^2)), 2)) %>% 
    dplyr::filter(var_exp != 0) %>% 
    dplyr::as_tibble()
  
  scree_plot <- ggplot2::ggplot(eigenvalues, ggplot2::aes(x = reorder(comp, -var_exp), y = var_exp)) +
    ggplot2::geom_col(ggplot2::aes(fill = var_exp), show.legend = FALSE) +
    ggplot2::labs(x = "Principal Component",
                  y = "% Variance Explained",
                  ) +
    theme_poma(axis_x_rotate = TRUE) +
    scale_fill_poma_c()
  
  # loading plot
  loadings_tb <- pca_res$rotation %>% 
    dplyr::as_tibble()
  
  loadings_plot <- loadings_tb %>%
    dplyr::mutate(feature = rownames(SummarizedExperiment::assay(data))) %>% 
    dplyr::arrange(dplyr::desc(abs(PC1))) %>%
    dplyr::select(feature, PC1:PC4) %>% 
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
  
  # scores
  score_data <- pca_res$x %>%
    dplyr::as_tibble()
  
  # biplot
  lam <- (pca_res$sdev[1:2] * sqrt(nrow(pca_res_df)))^load_length
  len <- t(t(pca_res$rotation[, 1:2]) * lam)*0.8
  pca_loadings <- data.frame(pca_res$rotation, to_x = len[,1], to_y = len[,2])
  
  biplot <- ggplot2::ggplot(pca_res_df, ggplot2::aes(x = PC1, y = PC2))+
    {if(ncol(SummarizedExperiment::colData(data)) != 0)ggplot2::geom_point(ggplot2::aes(fill = group), pch = 21, size = 3, alpha = 0.6)} +
    {if(ncol(SummarizedExperiment::colData(data)) == 0)ggplot2::geom_point(pch = 21, size = 3, alpha = 0.6)} +
    ggplot2::labs(x = paste0("PC1 (", round(100*(((pca_res$sdev^2)[1]) / sum(pca_res$sdev^2)), 2), "%)"),
                  y = paste0("PC2 (", round(100*(((pca_res$sdev^2)[2]) / sum(pca_res$sdev^2)), 2), "%)"),
                  fill = NULL,
                  color = NULL) +
    {if(ellipse)ggplot2::stat_ellipse(ggplot2::aes(color = group), type = "norm")} +
    {if(labels)ggplot2::geom_text(ggplot2::aes(label = sample_id), show.legend = FALSE)} +
    ggplot2::geom_segment(data = pca_loadings,
                          ggplot2::aes(x = 0, y = 0, xend = to_x, yend = to_y),
                          arrow = ggplot2::arrow(length = ggplot2::unit(1/2, "picas")), color = "grey19") +
    ggplot2::annotate("text", 
                      x = pca_loadings$to_x,
                      y = pca_loadings$to_y,
                      label = rownames(pca_loadings), size = 4) +
    theme_poma() +
    scale_fill_poma_d() +
    scale_color_poma_d()
  
  return(list(scree_plot = scree_plot, 
              scores_plot = scores_plot,
              scores = score_data, 
              eigenvalues = eigenvalues, 
              loadings = loadings_tb,
              loadings_plot = loadings_plot,
              biplot = biplot)
         )
}
    
