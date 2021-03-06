
#' Classical Heatmap
#'
#' @description This function returns a basic heatmap plot made with base R.
#' 
#' @param data A MSnSet object. First `pData` column must be the subject group/type.
#' @param sample_names Logical indicating if sample names should be plotted or not. Default is TRUE.
#' @param feature_names Logical indicating if feature names should be plotted or not. Default is FALSE.
#' @param show_legend Logical indicating if legend should be plotted or not. Default is TRUE.
#'
#' @export
#'
#' @return A heatmap.
#' @author Pol Castellano-Escuder
#'
#' @importFrom crayon red
#' @importFrom clisymbols symbol
#' @importFrom MSnbase pData exprs
#' @importFrom ComplexHeatmap HeatmapAnnotation Heatmap
#' 
#' @examples 
#' data("st000284")
#' 
#' st000284 %>% 
#'   PomaNorm() %>% 
#'   PomaHeatmap()
PomaHeatmap <- function(data, 
                        sample_names = TRUE,
                        feature_names = FALSE,
                        show_legend = TRUE){
  
  if (missing(data)) {
    stop(crayon::red(clisymbols::symbol$cross, "data argument is empty!"))
  }
  if(!is(data[1], "MSnSet")){
    stop(paste0(crayon::red(clisymbols::symbol$cross, "data is not a MSnSet object."), 
                " \nSee POMA::PomaMSnSetClass or MSnbase::MSnSet"))
  }
  
  total <- MSnbase::exprs(data)
  target <- MSnbase::pData(data)

  ha <- ComplexHeatmap::HeatmapAnnotation(df = data.frame(Group = target[,1]),
                                          show_legend = show_legend)

  ComplexHeatmap::Heatmap(total, name = "Value", 
                          top_annotation = ha,
                          show_row_names = feature_names, 
                          show_column_names = sample_names,
                          show_heatmap_legend = show_legend)
  
}

