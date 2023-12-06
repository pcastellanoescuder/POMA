.onAttach <- function(...) { # nocov start
  
  suppressWarnings({
    
    poma_welcome <- "Welcome to POMA!"
    ver <- paste("Version", utils::packageVersion("POMA"))
    poma_shiny <- "POMAShiny app: https://github.com/pcastellanoescuder/POMAShiny"
    
    packageStartupMessage(paste(poma_welcome, ver, poma_shiny, sep = "\n"))
    
  })
  
} # nocov end