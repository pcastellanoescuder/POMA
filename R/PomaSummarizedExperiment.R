
#' Create a `SummarizedExperiment` Object
#'
#' @description `PomaSummarizedExperiment` creates a `SummarizedExperiment` object from data frames.
#' 
#' @param metadata Metadata variables structured in columns. Sample ID must be the first column.
#' @param features Matrix of features. Each feature is a column.
#'
#' @export
#'
#' @return A `SummarizedExperiment` object.
#' @references Morgan M, Obenchain V, Hester J, Pagès H (2021). SummarizedExperiment: SummarizedExperiment container. R package version 1.24.0, https://bioconductor.org/packages/SummarizedExperiment.
#' @author Pol Castellano-Escuder
#'
#' @importFrom magrittr %>%
#' 
#' @examples 
#' data(iris)
#' 
#' # Create metadata: Data frame with sample names and a group factor
#' metadata <- data.frame(ID = 1:150, Group = iris$Species)
#' 
#' # Create features: `p` column data frame with features
#' features <- iris[, 1:4]
#' 
#' # Create a `SummarizedExperiment` object with `POMA`
#' object <- PomaSummarizedExperiment(metadata = metadata, features = features)
PomaSummarizedExperiment <- function(metadata = NULL,
                                     features = NULL,
                                     ...){
  
  if(missing(features)){
    stop("No features file")
  }
  
  # features
  features <- features %>% 
    as.data.frame() %>% 
    janitor::clean_names() %>% 
    dplyr::mutate_all(~ as.numeric(as.character(.)))
  
  # metadata
  if(!is.null(metadata)){
    if(!is.data.frame(metadata)){
      stop("metadata file is not a data.frame()")
    }
    if(sum(sapply(metadata, function(x)sum(is.na(x)))) > 0){
      stop("Missing values not allowed in metadata file")
    }
    if(nrow(metadata) != nrow(features)){
      stop("Different number of samples in metadata and features")
    }
    
    metadata <- metadata %>%
      as.data.frame() %>% 
      janitor::clean_names() %>% 
      tibble::remove_rownames() %>%
      dplyr::rename(sample_id = 1) %>%
      tibble::column_to_rownames("sample_id") %>% 
      dplyr::mutate_if(is.character, as.factor)

    rownames(features) <- rownames(metadata)
    
    # create a SummarizedExperiment object
    data <- SummarizedExperiment::SummarizedExperiment(assays = t(features), colData = metadata)
  } else {
    data <- SummarizedExperiment::SummarizedExperiment(assays = t(features))
  }
  
  if(validObject(data))
    return(data)
}

