
#' Colorectal Cancer Detection Using Targeted Serum Metabolic Profiling
#'
#' Colorectal cancer (CRC) is one of the most prevalent and deadly cancers in the world. Despite an expanding knowledge of its molecular 
#' pathogenesis during the past two decades, robust biomarkers to enable screening, surveillance, and therapy monitoring of CRC are still lacking. 
#' In this study, we present a targeted liquid chromatography-tandem mass spectrometry-based metabolic profiling approach for identifying biomarker 
#' candidates that could enable highly sensitive and specific CRC detection using human serum samples. In this targeted approach, 158 metabolites 
#' from 25 metabolic pathways of potential significance were monitored in 234 serum samples from three groups of patients (66 CRC patients, 
#' 76 polyp patients, and 92 healthy controls). Partial least squares-discriminant analysis (PLS-DA) models were established, which proved to 
#' be powerful for distinguishing CRC patients from both healthy controls and polyp patients. Receiver operating characteristic curves generated 
#' based on these PLS-DA models showed high sensitivities (0.96 and 0.89, respectively, for differentiating CRC patients from healthy controls 
#' or polyp patients); good specificities (0.80 and 0.88), and excellent areas under the curve (0.93 and 0.95) were also obtained. Monte Carlo 
#' cross validation (MCCV) was also applied, demonstrating the robust diagnostic power of this metabolic profiling approach.
#'
#' @format A `SummarizedExperiment` object: 224 samples, 113 metabolites, 4 covariables and 3 groups (CRC, Healthy and Polyp).
#' \describe{
#'   \item{metabolites}{113 serum metabolites.}
#'   \item{covariables}{Age at consent, Gender, Smoking Condition and Alcohol Consumption.}
#' }
#' @references Colorectal Cancer Detection Using Targeted Serum Metabolic Profiling, J. Proteome. Res., 2014, 13, 4120-4130.
#' @source \url{https://www.metabolomicsworkbench.org/data/DRCCMetadata.php?Mode=Study&StudyID=ST000284&StudyType=MS&ResultType=1%20target=_blank}
"st000284"


#' Targeted LC/MS of urine from boys with DMD and controls
#'
#' Duchenne Muscular Dystrophy (DMD) is an X-linked recessive form of muscular dystrophy that affects males via a mutation in the 
#' gene for the muscle protein, dystrophin. Progression of the disease results in severe muscle loss, ultimately leading to paralysis 
#' and death. Steroid therapy has been a commonly employed method for reducing the severity of symptoms. This study aims to quantify 
#' the urine levels of amino acids and organic acids in patients with DMD both with and without steroid treatment. Track the progression 
#' of DMD in patients who have provided multiple urine samples.
#'
#' @format A `SummarizedExperiment` object: 57 samples, 31 metabolites, 1 covariable and 2 groups (Controls and DMD).
#' \describe{
#'   \item{metabolites}{31 urine metabolites.}
#'   \item{covariables}{Steroid status.}
#' }
#' @source \url{https://www.metabolomicsworkbench.org/data/DRCCMetadata.php?Mode=Study&DataMode=AllData&StudyID=ST000336&StudyType=MS&ResultType=1#DataTabs}
"st000336"

