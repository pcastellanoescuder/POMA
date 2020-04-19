
#' Convert data frames to a simplified MSnSet Object
#'
#' @param target Metadata variables structured in columns. Sample ID must be the first column and group/type/treatment of the study must be the second column.
#' @param features Table of features. Each feature in one column.
#'
#' @export
#'
#' @return A simplified MSnSet object.
#' @references Laurent Gatto and Kathryn S. Lilley. MSnbase - an R/Bioconductor package for isobaric tagged mass spectrometry data visualization, processing and quantitation. Bioinformatics 28, 288-289 (2012).
#' @author Pol Castellano-Escuder
#'
#' @importFrom janitor clean_names
#' @importFrom tibble column_to_rownames remove_rownames
#' @importFrom magrittr %>%
#' @importFrom MSnbase MSnSet
PomaMSnSetClass <- function(target,
                            features){

  if(nrow(target) != nrow(features)){
    stop(crayon::red(clisymbols::symbol$cross, "Different number of samples between target and features!"))
  }
  if(missing(target)){
    stop(crayon::red(clisymbols::symbol$cross, "target required!"))
  }
  if(missing(features)){
    stop(crayon::red(clisymbols::symbol$cross, "features required!"))
  }
  if(!is.data.frame(target)){
    stop(crayon::red(clisymbols::symbol$cross, "target file is not a data.frame"))
  }

  source("R/simpleMSnSetClass.R")
  
  target <- as.data.frame(target)
  target <- remove_rownames(target)
  target <- target %>% janitor::clean_names()
  colnames(target)[1] <- "ID"
  target <- target %>% column_to_rownames("ID") %>% as.data.frame()

  features <- features %>% as.data.frame() %>% janitor::clean_names() %>% as.matrix()
  rownames(features) <- rownames(target)

  ## create a simpleMSnSet object

  # target <- new("AnnotatedDataFrame", data = target)
  # features <- new("ExpressionSet", exprs = t(features))
  # data <- new("simpleMSnSet", exprs = features, phenoData = target)
  data <- MSnbase::MSnSet(exprs = t(features), pData = target)

  if(validObject(data))
    return(data)

}

