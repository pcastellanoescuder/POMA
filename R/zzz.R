.onAttach <- function(...) { # nocov start
  
  suppressWarnings({
    
    orange <- crayon::make_style("orange", bg = FALSE)
    
    poma_welcome <- crayon::bold(orange("Welcome to POMA!"))
    ver <- paste(crayon::blue(clisymbols::symbol$tick), crayon::bold(orange(paste(utils::packageVersion("POMA"), "version"))))
    shiny1 <- paste(crayon::blue(clisymbols::symbol$tick), paste(crayon::bold(orange("POMA Shiny version: https://github.com/pcastellanoescuder/POMAShiny"))))
    # shiny2 <- paste(crayon::blue(clisymbols::symbol$tick), paste(crayon::bold(orange("POMAcounts Shiny module (for mass spectrometry spectral counts): http://uebshiny.vhir.org:3838/POMAcounts"))))
    info <- paste(crayon::blue(clisymbols::symbol$info), crayon::bold(paste(orange("For more detailed package information please visit https://pcastellanoescuder.github.io/POMA/"))))
    
    packageStartupMessage(paste(poma_welcome, ver, shiny1, shiny2, info, sep = "\n"))
    
  })
  
} # nocov end