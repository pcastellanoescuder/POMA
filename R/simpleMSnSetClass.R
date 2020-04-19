
## This file is a modification of the https://github.com/lgatto/MSnbase/blob/master/R/DataClasses.R file by Laurent Gatto.
## The aim is to create a simplified MSnSet class compatible with POMA and avoid all MSnbase package dependencies.
## Main reason is that "mzR" package difficults a lot the "POMA_Shiny" (https://github.com/pcastellanoescuder/POMA_Shiny) deployment on shinyapps.io

######################################################################
## simpleMSnProcess: Container for simpleMSnSet processing information
setClass("simpleMSnProcess",
         representation = representation(
           processing = "character",
           cleaned = "logical",
           normalised = "logical",
           MSnbaseVersion = "character"),
         contains = c("Versioned"),
         prototype  =  prototype(
           new("Versioned", versions = c(MSnProcess = "0.1.3")),
           processing = character(),
           MSnbaseVersion = character())) ## set in initialize()

#####################################################################
## The "simpleMSnSet" Class for MS Proteomics/Metabolomics Expression Data and Meta-Data
setClass("simpleMSnSet",
         representation = representation(
           processingData = "simpleMSnProcess"
           ),
         contains = "eSet",
         prototype = prototype(
           new("VersionedBiobase",
               versions = c(Biobase::classVersion("eSet"),
                            Biobase::classVersion("pSet"),
                            MSnSet = "0.4.0"))))

