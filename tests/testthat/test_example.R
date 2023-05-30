# run examples from nnSVG() function documentation
example(nnSVG, echo = FALSE)


test_that("example objects have correct class", {
  expect_s4_class(spe1, "SpatialExperiment")
  expect_s4_class(spe2, "SpatialExperiment")
})

test_that("example objects have correct dimensions", {
  expect_equal(dim(spe1), c(4, 3639))
  expect_equal(dim(spe2), c(4, 15003))
})

test_that("examples give correct numbers of significant SVGs", {
  expect_equal(as.numeric(table(rowData(spe1)$padj <= 0.05)), c(2, 2))
  expect_equal(as.numeric(table(rowData(spe2)$padj <= 0.05)), c(2, 2))
})

test_that("examples give correct p-values", {
  expect_equal(signif(rowData(spe1)$pval, 6), 
               c(0, 0, 1, 7.82181e-02))
  expect_equal(signif(rowData(spe2)$pval, 6), 
               c(2.25708e-13, 0, 4.65894e-01, 5.74463e-01))
})

