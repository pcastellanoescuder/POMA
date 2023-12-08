
test_that("PomaUMAP handles valid SummarizedExperiment objects", {
  data <- create_mock_summarized_experiment()
  umap_results <- PomaUMAP(data)
  expect_is(umap_results, "list")
  expect_true("umap_embeddings" %in% names(umap_results))
  expect_true("umap_plot" %in% names(umap_results))
})

test_that("PomaUMAP stops with non-SummarizedExperiment objects", {
  data <- data.frame(matrix(runif(100), ncol = 10))
  expect_error(PomaUMAP(data), "data is not a SummarizedExperiment object")
})

test_that("PomaUMAP works with different n_neighbors values", {
  data <- create_mock_summarized_experiment()
  umap_results <- PomaUMAP(data, n_neighbors = 5)
  expect_is(umap_results, "list")
})

test_that("PomaUMAP works with different n_components values", {
  data <- create_mock_summarized_experiment()
  umap_results <- PomaUMAP(data, n_components = 3)
  expect_is(umap_results, "list")
})

test_that("PomaUMAP provides expected output structure", {
  data <- create_mock_summarized_experiment()
  umap_results <- PomaUMAP(data)
  expect_is(umap_results$umap_embeddings, "tbl_df")
  expect_true(all(c("sample", "UMAP1", "UMAP2", "clust", "member_prob") %in% names(umap_results$umap_embeddings)))
})

