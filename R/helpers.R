
help_extract <- function(fun, 
                         section = "Description",
                         ...) {
  x <- capture.output(tools:::Rd2txt(utils:::.getHelpFile(help(fun, ...)), options = list(sectionIndent = 0)))
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

title_extract <- function(fun, 
                          ...) {
  x <- capture.output(tools:::Rd2txt(utils:::.getHelpFile(help(fun, ...)), options = list(sectionIndent = 0)))
  x <- gsub("_\b", "", x, fixed = TRUE)
  title <- x[1]
  return(title)
}

make_legend <- function(fun, 
                        html = FALSE, 
                        ...) {
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

