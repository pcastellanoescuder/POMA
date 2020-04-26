## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  dpi = 200,
  fig.align = "center"
)

## ---- eval = FALSE------------------------------------------------------------
#  # install.packages("devtools")
#  devtools::install_github("pcastellanoescuder/POMA")

## ---- warning = FALSE, message = FALSE, comment = FALSE-----------------------
library(POMA)
library(ggplot2) # to tune POMA default plots
library(patchwork) # to show plots together
library(tidyverse) # it's just tidyverse;)

## -----------------------------------------------------------------------------
# load example data
data("st000336")

# # subset example data
# my_features <- t(exprs(st000284))[c(1:8, 120:127) ,]
# 
# my_target <- pData(st000284)[c(1:8, 120:127) ,]
# my_target <- my_target %>% rownames_to_column("ID")
# 
# # create a smaller MSnSet object from example data
# example_data <- PomaMSnSetClass(features = my_features, target = my_target)

## -----------------------------------------------------------------------------
imputed <- PomaImpute(st000336, ZerosAsNA = T, RemoveNA = T, cutoff = 20, method = "knn")
imputed

## -----------------------------------------------------------------------------
normalized <- PomaNorm(imputed, method = "log_pareto")
normalized

## ---- fig.height = 4, fig.width = 8, message = FALSE, comment = FALSE---------
p1 <- PomaBoxplots(imputed, group = "samples", jitter = FALSE) +
  ggtitle("Not Normalized") +
  theme(legend.position = "none") # data before normalization

p2 <- PomaBoxplots(normalized, group = "samples", jitter = FALSE) +
  ggtitle("Normalized") # data after normalization

p1 + p2

## ---- fig.height = 4, fig.width = 8, message = FALSE, comment = FALSE---------
p3 <- PomaDensity(imputed, group = "features") +
  ggtitle("Not Normalized") +
  theme(legend.position = "none") # data before normalization

p4 <- PomaDensity(normalized, group = "features") +
  ggtitle("Normalized") # data after normalization

p3 + p4

## -----------------------------------------------------------------------------
ttest_res <- PomaUnivariate(normalized, method = "ttest", 
                            paired = FALSE, var_equal = FALSE, adjust = "fdr")
ttest_res %>% 
  rownames_to_column() %>% 
  head()

## -----------------------------------------------------------------------------
limma_res <- PomaLimma(normalized, contrast = "Controls-DMD", covariates = FALSE, adjust = "fdr")
limma_res %>% 
  rownames_to_column() %>% 
  head()

## -----------------------------------------------------------------------------
multiv_pca <- PomaMultivariate(normalized, method = "pca", components = 5,
                               scale = FALSE,
                               center = FALSE)

## -----------------------------------------------------------------------------
multiv_pca$score_data %>% head()

## ---- fig.height = 4, fig.width = 10------------------------------------------
p5 <- multiv_pca$screeplot +
  ggtitle("Scree Plot")

p6 <- multiv_pca$scoresplot +
  ggtitle("Scores Plot")

p5 + p6

## -----------------------------------------------------------------------------
rf_res <- PomaRandForest(normalized, folds = 3, nvar = 10)

## -----------------------------------------------------------------------------
rf_res$confusion_matrix

## ---- fig.height = 4, fig.width = 8-------------------------------------------
rf_res$gini_plot

