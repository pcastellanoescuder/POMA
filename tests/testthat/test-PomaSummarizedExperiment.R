context("PomaSummarizedExperiment")

test_that("PomaSummarizedExperiment works", {

  target <- data.frame(ID = c("One", "Two", "Three", "Four"), Group = c("Trtd", "Ctrl", "Trtd", "Ctrl"), Smoking = c(1,0,0,1))
  target2 <- data.frame(ID = c("Five", "One", "Three", "Two"), Group = c("Ctrl", "Trtd", "Trtd", "Ctrl"), Smoking = c(0,0,0,1))
  target_error <- as.matrix(target)
  target_error_2 <- data.frame(ID = c("Five", "One", "Three", "Two"), Group = c("Ctrl", "Trtd", "Trtd", "Ctrl"), Smoking = c(0,0,NA,1))
    
  features <- data.frame(Feat.1 = c(1,2,3,4), Feat.2 = c(6,3,7,3), Feat.3 = c(3,5,23,24))
  features_error <- data.frame(Feat.1 = c(1,2,3,4,5), Feat.2 = c(6,3,7,4,3), Feat.3 = c(3,4,5,23,24))

  a <- PomaSummarizedExperiment(target, features)
  b <- PomaSummarizedExperiment(target2, features)

  ##

  expect_true(validObject(a))
  expect_true(validObject(b))

  ##

  expect_false(all(rownames(SummarizedExperiment::colData(a)) == rownames(SummarizedExperiment::colData(b))))

  ##

  expect_error(PomaSummarizedExperiment(target_error, features))
  expect_error(PomaSummarizedExperiment(target, features_error))
  expect_error(PomaSummarizedExperiment(target))
  expect_error(PomaSummarizedExperiment(features))
  
  ##
  
  expect_equal(colnames(features), rownames(SummarizedExperiment::assay(a)))
  expect_false(all(colnames(target)[2:3] == SummarizedExperiment::colData(a)))
  expect_false(all(colnames(target2)[2:3] == SummarizedExperiment::colData(b)))

  ##
  
  expect_error(PomaSummarizedExperiment(target_error_2, features))
  
})

