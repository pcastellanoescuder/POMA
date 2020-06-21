
#' Classical Heatmap
#'
#' @description This function returns a basic heatmap plot made with base R.
#' 
#' @param data A MSnSet object. First `pData` column must be the subject group/type.
#' @param sample_names Logical indicating if sample names should be plotted or not. Default is TRUE.
#' @param feature_names Logical indicating if feature names should be plotted or not. Default is FALSE.
#'
#' @export
#'
#' @return A heatmap.
#' @author Pol Castellano-Escuder
#'
#' @importFrom crayon red
#' @importFrom clisymbols symbol
#' @importFrom Biobase pData exprs
#' @importFrom ComplexHeatmap HeatmapAnnotation Heatmap
PomaHeatmap <- function(data, 
                        sample_names = TRUE,
                        feature_names = FALSE){
  
  if (missing(data)) {
    stop(crayon::red(clisymbols::symbol$cross, "data argument is empty!"))
  }
  if(!is(data)[1] == "MSnSet"){
    stop(paste0(crayon::red(clisymbols::symbol$cross, "data is not a MSnSet object."), 
                " \nSee POMA::PomaMSnSetClass or MSnbase::MSnSet"))
  }
  
  total <- Biobase::exprs(data)
  target <- Biobase::pData(data)

  ha <- ComplexHeatmap::HeatmapAnnotation(df = data.frame(Group = target[,1]))

  ComplexHeatmap::Heatmap(total, name = "Value", top_annotation = ha,
                          show_row_names = feature_names, show_column_names = sample_names)
  
}

