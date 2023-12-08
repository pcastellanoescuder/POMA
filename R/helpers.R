
help_extract <- function(fun, 
                         section = "Description") {
  x <- capture.output(tools::Rd2txt(utils:::.getHelpFile(help(fun, ...)), options = list(sectionIndent = 0)))
  B <- grep("^_", x)
  x <- gsub("_\b", "", x, fixed = TRUE)
  X <- rep(FALSE, length(x))
  X[B] <- 1
  out <- split(x, cumsum(X))
  out <- out[[which(sapply(out, function(x) grepl(section, x[1], fixed = TRUE)))]][-c(1, 2)]
  while (TRUE) {
    out <- out[-length(out)]
    if (out[length(out)] != "") {
      break
    }
  }
  if (section == "Description") {
    out <- paste0(out, collapse = " ")
  }
  return(out)
}

title_extract <- function(fun) {
  x <- capture.output(tools::Rd2txt(utils:::.getHelpFile(help(fun, ...)), options = list(sectionIndent = 0)))
  x <- gsub("_\b", "", x, fixed = TRUE)
  title <- x[1]
  return(title)
}

make_legend <- function(fun, 
                        html = FALSE) {
  description <- help_extract(fun, package = POMA, section = "Description")
  title <- title_extract(fun, package = POMA)
  if (html) {
    legend <- paste0("<b>", title, ".</b> ", description)
  }
  else {
    legend <- paste0(title, ". ", description)
  }
  return(legend)
}

create_mock_summarized_experiment <- function(binary = FALSE, paired = FALSE, integers = FALSE) {
  
  if (!binary) {
    g_labels <- sample(c("A", "B", "C"), 20, replace = TRUE)
  } else {
    g_labels <- sample(c("A", "B"), 20, replace = TRUE)
    if (paired) {
      g_labels <- c(rep("A", 10), rep("B", 10))
    }
  }
  
  if (integers) {
    matrix_data <- matrix(sample(1:100, 20 * 10, replace = TRUE), nrow = 20, ncol = 10)
  } else {
    matrix_data <- matrix(runif(100), nrow = 20)
  }
  
  col_data <- data.frame(sample = paste0("Sample", 1:20), group = g_labels)
  PomaCreateObject(features = matrix_data, metadata = col_data)
}

create_mock_data <- function() {
  features <- as.data.frame(matrix(runif(100), ncol = 10))
  metadata <- data.frame(ID = 1:10, Group = factor(rep(c("A", "B"), each = 5)))
  list(features = features, metadata = metadata)
}

