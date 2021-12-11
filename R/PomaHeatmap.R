
#' Classical Heatmap
#'
#' @description This function returns a basic heatmap plot made with base R.
#' 
#' @param data A SummarizedExperiment object. First `colData` column must be the subject group/type.
#' @param cols Numerical vector indicating the column index of variables in `colData` to be displayed. Default is 1 (main group).
#' @param sample_names Logical indicating if sample names should be plotted or not. Default is TRUE.
#' @param feature_names Logical indicating if feature names should be plotted or not. Default is FALSE.
#' @param show_legend Logical indicating if legend should be plotted or not. Default is TRUE.
#'
#' @export
#'
#' @return A heatmap.
#' @author Pol Castellano-Escuder
#'
#' @importFrom SummarizedExperiment assay colData
#' @importFrom ComplexHeatmap HeatmapAnnotation Heatmap
#' 
#' @examples 
#' data("st000284")
#' 
#' st000284 %>% 
#'   PomaNorm() %>% 
#'   PomaHeatmap()
PomaHeatmap <- function(data, 
                        cols = 1,
                        sample_names = TRUE,
                        feature_names = FALSE,
                        show_legend = TRUE){
  
  if (missing(data)) {
    stop("data argument is empty!")
  }
  if(!is(data[1], "SummarizedExperiment")){
    stop("data is not a SummarizedExperiment object. \nSee POMA::PomaSummarizedExperiment or SummarizedExperiment::SummarizedExperiment")
  }
  
  total <- SummarizedExperiment::assay(data)
  target <- SummarizedExperiment::colData(data)

  df <- data.frame(target[, cols])
  colnames(df) <- names(target[cols])
  ha <- ComplexHeatmap::HeatmapAnnotation(df = df,
                                          show_legend = show_legend)

  suppressMessages(
    ComplexHeatmap::Heatmap(total, 
                            name = "Value", 
                            top_annotation = ha,
                            show_row_names = feature_names, 
                            show_column_names = sample_names,
                            show_heatmap_legend = show_legend)
  )
  
}

