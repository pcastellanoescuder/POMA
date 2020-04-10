context("PomaMSnSetClass")

test_that("PomaMSnSetClass works", {

  target <- data.frame(ID = c("One", "Two", "Three", "Four"), Group = c("Trtd", "Ctrl", "Trtd", "Ctrl"), Smoking = c(1,0,0,1))
  target2 <- data.frame(ID = c("Five", "One", "Three", "Two"), Group = c("Ctrl", "Trtd", "Trtd", "Ctrl"), Smoking = c(0,0,0,1))
  target_error <- as.matrix(target)

  features <- data.frame(feat1 = c(1,2,3,4), feat2 = c(6,3,7,3), feat3 = c(3,5,23,24))
  features_error <- data.frame(feat1 = c(1,2,3,4,5), feat2 = c(6,3,7,4,3), feat3 = c(3,4,5,23,24))

  a <- PomaMSnSetClass(target, features)
  b <- PomaMSnSetClass(target2, features)

  ##

  expect_true(validObject(a))
  expect_true(validObject(b))

  ##

  expect_false(all(Biobase::sampleNames(a) == Biobase::sampleNames(b)))

  ##

  expect_error(PomaMSnSetClass(target_error, features))
  expect_error(PomaMSnSetClass(target, features_error))
  expect_error(PomaMSnSetClass(target))
  expect_error(PomaMSnSetClass(features))

})
