# POMA 1.5.17

* `MSnbase::MSnSet` class has been replaced by the `SummarizedExperiment` class
* Color scale for all plots set to `viridis` (without yellow)
* All output tables provided as `tibble` istead of `matrix` or `data.frame`
* Allow users to select specific covariates and their position (importance) in the model for `PomaLimma`, `PomaUnivariate(method = "anova")`, and `PomaOddsRatio`
* Add SD to `PomaUnivariate` output tables
* Reduce dependencies
* Update vignettes
* Update documentation
* Some other major and minor improvements
* Some minor bugs and typos fixed
* New "loading plot" in PCA
* Compute p-values and FDR in `PomaCorr`

# POMA 1.1.15

* Remove `reshape2` and `Biobase` packages from Imports
* Implement viridis palette for `PomaBoxplots`, `PomaDensity`, and `PomaMultivariate`
* Update `mixOmics` output names in `PomaMultivariate`
* New package description and biocViews
* Bug fixed in `PomaMultivariate`

# POMA 1.1.8

* `PomaNorm` and `PomaImpute` warnings when methos parameter is missing passed to messages
* Minor bugs fixed
* Minor changes in plots style
* `plotly` package used in `PomaVolcano` switched to Suggests
* _Bioconductor_ badge table added to README

# POMA 1.0.0

* Released to Bioconductor 3.12

# POMA 0.99.45

* PomaOutliers bug fixed
* New references in PomaRankProd help
* PomaLimma and PomaUnivariate bug related with one covariate analysis fixed
* The elbow method to calculate the optimum number of clusters has been added in PomaClust function

# POMA 0.99.37

* POMA has been accepted to Bioconductor!
* PomaRankProd minor bug fixed
* Authors updated
* Bioconductor logo

# POMA 0.99.33

* All Bioconductor issues in the review process have been addressed
* POMA EDA vignette added 

# POMA 0.99.16

* BiocCeck requirements fixed
* Examples added in functions
* pkgdown files removed from master branch

# POMA 0.99.0

* POMA is now submitted to Bioconductor!
* All features implemented work as expected
* All tests finished
* Achieved desired coverage (>95%)
* Vignettes and documentation are ready for the first release

