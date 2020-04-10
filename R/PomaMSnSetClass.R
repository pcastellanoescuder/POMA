
#' Convert data frames to an MSnSet Object
#'
#' @param target Metadata variables structured in columns. Sample ID must be the first column and group/type/treatment of the study must be the second column.
#' @param features Table of features. Each feature in one column.
#'
#' @export
#'
#' @return An MSnSet object.
#' @references Laurent Gatto and Kathryn S. Lilley. MSnbase - an R/Bioconductor package for isobaric tagged mass spectrometry data visualization, processing and quantitation. Bioinformatics 28, 288-289 (2012).
#' @author Pol Castellano-Escuder
#'
#' @importFrom MSnbase MSnSet
#' @importFrom tibble column_to_rownames remove_rownames
#' @importFrom magrittr %>%
PomaMSnSetClass <- function(target,
                            features){

  if(nrow(target) != nrow(features)){
    stop(crayon::red(clisymbols::symbol$cross, "Different number of sambles between target and features!"))
  }
  if(missing(target)){
    stop(crayon::red(clisymbols::symbol$cross, "target required!"))
  }
  if(missing(features)){
    stop(crayon::red(clisymbols::symbol$cross, "features required!"))
  }
  if(!is.data.frame(target)){
    stop(crayon::red(clisymbols::symbol$cross, "target must be a data frame!"))
  }

  target <- as.data.frame(target)
  target <- remove_rownames(target)
  colnames(target)[1] <- "ID"

  target <- target %>% column_to_rownames("ID")

  features <- as.matrix(features)
  rownames(features) <- rownames(target)

  data <- MSnbase::MSnSet(exprs = t(features), pData = target)

  if(validObject(data))
    return(data)

}

