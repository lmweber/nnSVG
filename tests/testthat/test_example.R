# run example from nnSVG function documentation
example(nnSVG, echo = FALSE)

test_that("example has correct class", {
  expect_s4_class(spe, "SpatialExperiment")
})

test_that("example has correct dimensions", {
  expect_equal(dim(spe), c(4, 3639))
})

test_that("example identifies correct number of significant SVGs", {
  expect_equal(as.numeric(table(rowData(spe)$padj <= 0.05)), c(1, 3))
})

