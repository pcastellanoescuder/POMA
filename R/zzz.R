.onAttach <- function(...) { # nocov start
  
  suppressWarnings({
    
    poma_welcome <- "Welcome to POMA!"
    ver <- paste("Version", utils::packageVersion("POMA"))
    poma_shiny <- "POMAShiny app: https://github.com/pcastellanoescuder/POMAShiny"
    info <- "For more detailed package information please visit https://pcastellanoescuder.github.io/POMA/"
    
    packageStartupMessage(paste(poma_welcome, ver, poma_shiny, info, sep = "\n"))
    
  })
  
} # nocov end