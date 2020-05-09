
#' Classical Heatmap
#'
#' @description This function returns a basic heatmap plot made with base R.
#' 
#' @param data A MSnSet object. First `pData` column must be the subject group/type.
#' @param scale Logical that indicates if heatmap values should be scaled. Default is TRUE.
#' @param scale_by Variable to scale heatmap values. Options are "samples" and "features". Option "samples" (default) will `scale()` each sample considering all features and option "features" will `scale()` each feature considering all samples.
#'
#' @export
#'
#' @return A heatmap.
#' @author Pol Castellano-Escuder
#'
#' @importFrom crayon red
#' @importFrom clisymbols symbol
#' @importFrom Biobase pData exprs
PomaHeatmap <- function(data, 
                        scale = TRUE,
                        scale_by = "features"){
  
  if (missing(data)) {
    stop(crayon::red(clisymbols::symbol$cross, "data argument is empty!"))
  }
  if(!(class(data) == "MSnSet")){
    stop(paste0(crayon::red(clisymbols::symbol$cross, "data is not a MSnSet object."), 
                " \nSee POMA::PomaMSnSetClass or MSnbase::MSnSet"))
  }
  if (!(scale_by %in% c("samples", "features"))) {
    stop(crayon::red(clisymbols::symbol$cross, "Incorrect value for scale_by argument!"))
  }
  
  total <- t(Biobase::exprs(data))
  target <- Biobase::pData(data)

  my_group <- as.numeric(as.factor(target[,1]))
  colSide <- topo.colors(length(table(my_group)))[my_group]
  colMain <- colorRampPalette( c("green", "black", "red"), space = "rgb")(100)

  if(scale){
    if(scale_by == "samples"){
      heatmap(scale(t(total)), ColSideColors = colSide, col = colMain)
    }
    if(scale_by == "features"){
      heatmap(t(scale(total)), ColSideColors = colSide, col = colMain)
    }
  }
  else{
    heatmap(t(total), ColSideColors = colSide, col = colMain)
  }
}

