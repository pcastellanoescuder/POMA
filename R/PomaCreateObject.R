
#' Detect decimals
#'
#' Detect decimal variables.
#' 
#' @param data A data matrix (samples in rows).
detect_decimals <- function(data) {
  decimal_columns <- sapply(data, function(col) {
    any(col %% 1 != 0)
  })
  return(decimal_columns)
}

#' Create a `SummarizedExperiment` Object
#'
#' @description `PomaCreateObject` creates a `SummarizedExperiment` object from data frames.
#' 
#' @param metadata Metadata variables structured in columns. Sample ID must be the first column.
#' @param features Matrix of features. Each feature is a column.
#' @param factor_levels Numeric. Integer variables with more levels than indicated by this parameter will be treated as factors.
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
#' object <- PomaCreateObject(metadata = metadata, features = features)
PomaCreateObject <- function(metadata = NULL,
                             features = NULL,
                             factor_levels = 10){
  
  if(missing(features)){
    stop("No features file")
  }
  
  # features
  features <- features %>% 
    as.data.frame() %>% 
    janitor::clean_names(case = "all_caps") %>% 
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

    numeric_vars <- metadata %>% 
      dplyr::select_if(is.numeric)
    
    if (ncol(numeric_vars) > 0) {
      decimal_vars <- detect_decimals(numeric_vars)
      decimal_vars <- names(decimal_vars)[decimal_vars]
      
      numerical_float <- numeric_vars %>% 
        dplyr::select(dplyr::all_of(decimal_vars)) %>% 
        dplyr::mutate_all(as.double)
      
      numerical_integer <- numeric_vars %>% 
        dplyr::select(-dplyr::all_of(decimal_vars)) %>% 
        dplyr::mutate_all(as.integer)
      
      if (any(apply(numerical_integer, 2, function(x) length(table(x))) > factor_levels)) {
        other_decimals <- names(which(apply(numerical_integer, 2, function(x) length(table(x))) > factor_levels))
        decimal_vars <- c(decimal_vars, other_decimals)
        
        numerical_float <- numeric_vars %>% 
          dplyr::select(dplyr::all_of(decimal_vars)) %>% 
          dplyr::mutate_all(as.double)
        
        numerical_integer <- numeric_vars %>% 
          dplyr::select(-dplyr::all_of(decimal_vars)) %>% 
          dplyr::mutate_all(as.integer)
      }
      
      metadata <- metadata %>% 
        dplyr::select(-dplyr::all_of(colnames(numerical_integer))) %>% 
        dplyr::select(-dplyr::all_of(colnames(numerical_float))) %>% 
        dplyr::bind_cols(numerical_integer) %>% 
        dplyr::bind_cols(numerical_float) %>% 
        dplyr::mutate_if(is.integer, as.factor) 
    }
    
    rownames(features) <- rownames(metadata)
    
    # create a SummarizedExperiment object
    data <- SummarizedExperiment::SummarizedExperiment(assays = t(features), colData = metadata)
  } else {
    data <- SummarizedExperiment::SummarizedExperiment(assays = t(features))
  }
  
  if(validObject(data))
    return(data)
}

