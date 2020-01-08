PomaMSnSetConstructor <- function(exprs, pData, fData){

  if (missing(exprs)) {
    stop("exprs argument is empty!")
  }
  if (missing(pData)) {
    stop("pData argument is empty!")
  }
  if(!missing(fData)){
    fData <- as.data.frame(fData)
  }

  metabolites <- t(exprs)
  pData <- as.data.frame(pData)
  pData <- tibble::column_to_rownames(pData, "ID")

  if (nrow(pData) != ncol(metabolites)) {
    stop("The number of rows in pData and exprs is not equal!")
  }

  m <- MSnbase::MSnSet(exprs = metabolites, pData = pData, fData = fData)

}

