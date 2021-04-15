context("PomaMSnSetClass")

test_that("PomaMSnSetClass works", {

  target <- data.frame(ID = c("One", "Two", "Three", "Four"), Group = c("Trtd", "Ctrl", "Trtd", "Ctrl"), Smoking = c(1,0,0,1))
  target2 <- data.frame(ID = c("Five", "One", "Three", "Two"), Group = c("Ctrl", "Trtd", "Trtd", "Ctrl"), Smoking = c(0,0,0,1))
  target_error <- as.matrix(target)
  target_error_2 <- data.frame(ID = c("Five", "One", "Three", "Two"), Group = c("Ctrl", "Trtd", "Trtd", "Ctrl"), Smoking = c(0,0,NA,1))
    
  features <- data.frame(Feat.1 = c(1,2,3,4), Feat.2 = c(6,3,7,3), Feat.3 = c(3,5,23,24))
  features_error <- data.frame(Feat.1 = c(1,2,3,4,5), Feat.2 = c(6,3,7,4,3), Feat.3 = c(3,4,5,23,24))

  a <- PomaMSnSetClass(target, features)
  b <- PomaMSnSetClass(target2, features)

  ##

  expect_true(validObject(a))
  expect_true(validObject(b))

  ##

  expect_false(all(MSnbase::sampleNames(a) == MSnbase::sampleNames(b)))

  ##

  expect_error(PomaMSnSetClass(target_error, features))
  expect_error(PomaMSnSetClass(target, features_error))
  expect_error(PomaMSnSetClass(target))
  expect_error(PomaMSnSetClass(features))
  
  ##
  
  expect_equal(colnames(features), colnames(t(MSnbase::exprs(a))))
  expect_false(all(colnames(target)[2:3] == colnames(MSnbase::pData(a))))
  expect_false(all(colnames(target2)[2:3] == colnames(MSnbase::pData(b))))

  ##
  
  expect_error(PomaMSnSetClass(target_error_2, features))
  
})

