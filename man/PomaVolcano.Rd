% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PomaVolcano.R
\name{PomaVolcano}
\alias{PomaVolcano}
\title{Volcano Plot}
\usage{
PomaVolcano(
  data,
  pval = "raw",
  pval_cutoff = 0.05,
  adjust = "fdr",
  log2FC = 0.6,
  xlim = 2,
  labels = FALSE,
  paired = FALSE,
  var_equal = FALSE,
  interactive = FALSE,
  plot_title = TRUE
)
}
\arguments{
\item{data}{A MSnSet object. First \code{pData} column must be the subject group/type. Only for two group data!}

\item{pval}{Select a pvalue type to generate the volcano plot. Options are: "raw" and "adjusted".}

\item{pval_cutoff}{Numeric. Define the pvalue cutoff (horizontal line).}

\item{adjust}{Multiple comparisons correction method for t test result. Options are: "fdr", "holm", "hochberg", "hommel", "bonferroni", "BH" and "BY".}

\item{log2FC}{Numeric. Define the log2 fold change cutoff (vertical lines).}

\item{xlim}{Numeric. Define the limits for x axis.}

\item{labels}{Logical that indicates if selected labels will be plotted or not. Defaul is FALSE.}

\item{paired}{Logical that indicates if the data is paired or not.}

\item{var_equal}{Logical that indicates if the data variance is equal or not.}

\item{interactive}{Logical that indicates if an interactive plot will be plotted or not. Defaul is FALSE.}

\item{plot_title}{Logical that indicates if title will be plotted or not. Defaul is TRUE.}
}
\value{
A ggplot2 object.
}
\description{
PomaVolcano() generates a volcano plot from the PomaUnivariate(method = "ttest") result. The data can't have negative values.
}
\examples{
data("st000336")

st000336 \%>\% 
  PomaImpute() \%>\%
  PomaVolcano()
}
\author{
Pol Castellano-Escuder
}
