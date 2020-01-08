context("PomaMSnSetConstructor")

test_that("PomaMSnSetConstructor works", {

  data <- readr::read_csv("data/MET_CRC_ST000284.csv")
  target <- readr::read_csv("data/TAR_CRC_ST000284.csv")

  m <- PomaMSnSetConstructor(exprs = data, pData =  target)

  # exprs(m)
  # pData(m)
  # fData(x)
  # sampleNames(m)
  # varLabels(m)
  # featureNames(m)

}


