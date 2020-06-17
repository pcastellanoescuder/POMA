FROM r-base:3.6.0

MAINTAINER Pol Castellano-Escuder <polcaes@gmail.com>

LABEL authors = "polcaes@gmail.com" \
      description = "Docker image containing the POMA R package latest version from GitHub"

# Install POMA dependencies

## CRAN

RUN R -e "install.packages(c('remotes', 'reshape2', 'ggplot2', 'tidyr', 'dplyr', 'tibble', 'stringr', 'caret', 'crayon', 'clisymbols', 'patchwork', 'knitr', 'rmarkdown', 'ggrepel', 'magrittr', 'randomForest', 'broom', 'glmnet', 'plotly', 'ggcorrplot', 'glasso', 'ggraph', 'vegan', 'BiocManager'), 
repos='http://cran.rstudio.com/')"

## BIOCONDUCTOR

RUN R -e "BiocManager::install(c('impute', 'RankProd', 'mixOmics', 'MSnbase', 'limma', 'Biobase', 'ComplexHeatmap'))"

# Install POMA from GitHub

RUN installGithub.r "pcastellanoescuder/POMA"

COPY POMA_*.tar.gz /app.tar.gz
RUN remotes::install_local('/app.tar.gz')
CMD R -e 'library(POMA)'
