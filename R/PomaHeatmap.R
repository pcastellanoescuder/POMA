
#' Heatmap Plot
#'
#' @description `PomaHeatmap` generates a heatmap. 
#' 
#' @param data A `SummarizedExperiment` object.
#' @param covs Character vector. Indicates the names of `colData` columns to be included as covariates. Default is NULL (no covariates).
#' @param sample_names Logical. Indicates if sample names should be displayed or not. Default is TRUE.
#' @param feature_names Logical. Indicates if feature names should be displayed or not. Default is FALSE.
#' @param show_legend Logical. Indicates if legend should be displayed or not. Default is TRUE.
#'
#' @export
#'
#' @return A heatmap plot.
#' @author Pol Castellano-Escuder
#' 
#' @examples 
#' data("st000284")
#' 
#' # Basic heatmap
#' st000284 %>% 
#'   PomaNorm() %>% 
#'   PomaHeatmap()
#'   
#' # Heatmap with one covariate  
#' st000284 %>% 
#'   PomaNorm() %>% 
#'   PomaHeatmap(covs = "factors")
#'   
#' # Heatmap with two covariates
#' st000284 %>% 
#'   PomaNorm() %>% 
#'   PomaHeatmap(covs = c("factors", "smoking_condition"))
PomaHeatmap <- function(data, 
                        covs = NULL,
                        sample_names = TRUE,
                        feature_names = FALSE,
                        show_legend = TRUE){
  
  if(!is(data, "SummarizedExperiment")){
    stop("data is not a SummarizedExperiment object. \nSee POMA::PomaCreateObject or SummarizedExperiment::SummarizedExperiment")
  }
  
  plot_data <- SummarizedExperiment::assay(data)
  
  if (ncol(SummarizedExperiment::colData(data)) > 0 & !is.null(covs)) {
    metadata <- SummarizedExperiment::colData(data) %>%
      as.data.frame() %>%
      dplyr::select(dplyr::all_of(covs))
    
    heatmap_annotations <- ComplexHeatmap::HeatmapAnnotation(df = metadata, show_legend = show_legend)
    
    suppressMessages(
      ComplexHeatmap::Heatmap(plot_data, 
                              name = "Value", 
                              top_annotation = heatmap_annotations,
                              show_row_names = feature_names, 
                              show_column_names = sample_names,
                              show_heatmap_legend = show_legend)
    )
  } else {
    suppressMessages(
      ComplexHeatmap::Heatmap(plot_data, 
                              name = "Value", 
                              show_row_names = feature_names, 
                              show_column_names = sample_names,
                              show_heatmap_legend = show_legend)
    )
  }
}

