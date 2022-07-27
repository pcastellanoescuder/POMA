
#' Convert data frames to a SummarizedExperiment Object
#'
#' @description This function converts data frame objects to a SummarizedExperiment object.
#' 
#' @param target Metadata variables structured in columns. Sample ID must be the first column and group/type/treatment of the study must be the second column.
#' @param features Table of features. Each feature in one column.
#'
#' @export
#'
#' @return A SummarizedExperiment object.
#' @references Morgan M, Obenchain V, Hester J, Pagès H (2021). SummarizedExperiment: SummarizedExperiment container. R package version 1.24.0, https://bioconductor.org/packages/SummarizedExperiment.
#' @author Pol Castellano-Escuder
#'
#' @importFrom magrittr %>%
#' 
#' @examples 
#' data(iris)
#' 
#' # create target: two column (or more) data frame with IDs and Group factor
#' target <- data.frame(ID = 1:150, Group = iris$Species)
#' 
#' # create features: p column data frame (or matrix) with features
#' features <- iris[,1:4]
#' 
#' # create an SummarizedExperiment object with POMA
#' object <- PomaSummarizedExperiment(target = target, features = features)
PomaSummarizedExperiment <- function(target,
                                     features){

  if(nrow(target) != nrow(features)){
    stop( "Different number of samples between target and features!")
  }
  if(missing(target)){
    stop("target required!")
  }
  if(missing(features)){
    stop("features required!")
  }
  if(!is.data.frame(target)){
    stop("target file is not a data.frame")
  }
  if(sum(sapply(target, function(x)sum(is.na(x)))) > 0){
    stop("missing values not allowed in target file")
  }
  
  target <- target %>%
    dplyr::as_tibble() %>% 
    tibble::remove_rownames() %>%
    dplyr::rename(ID = 1, group = 2) %>% 
    dplyr::mutate(group = as.factor(group)) %>%
    tibble::column_to_rownames("ID")

  features <- features %>% 
    as.data.frame()
  features <- as.matrix(sapply(features, function(x)as.numeric(as.character(x))))
  rownames(features) <- rownames(target)

  ## create a SummarizedExperiment object
  data <- SummarizedExperiment::SummarizedExperiment(assays = t(features), colData = target)
  data@NAMES <- make.names(data@NAMES)
  
  if(validObject(data))
    return(data)

}

