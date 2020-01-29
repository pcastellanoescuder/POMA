
#' Convert data frames to an MSnSet Object
#'
#' @param target Metadata variables structured in columns. Sample ID must be the first column and group/type of study must be the second column.
#' @param features Table of features. Each feature in one column.
#'
#' @export
#'
#' @return An MSnSet object.
#' @references Laurent Gatto and Kathryn S. Lilley. MSnbase - an R/Bioconductor package for isobaric tagged mass spectrometry data visualization, processing and quantitation. Bioinformatics 28, 288-289 (2012).
#' @author Pol Castellano-Escuder
#'
#' @importFrom MSnbase MSnSet
#' @importFrom tibble column_to_rownames
PomaMSnSetClass <- function(target, features){

  target <- as.data.frame(target)
  colnames(target)[1] <- "ID"
  target <- column_to_rownames(target, "ID")

  rownames(features) <- rownames(target)

  data <- MSnbase::MSnSet(exprs = t(features), pData = target)

  if (validObject(data))
    return(data)

}

