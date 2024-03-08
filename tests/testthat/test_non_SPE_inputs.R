# create random data
set.seed(123)
input <- matrix(runif(1000), nrow = 10)
spatial_coords <- matrix(runif(200), ncol = 2)

# run nnSVG with non-SpatialExperiment inputs
set.seed(123)
out <- nnSVG(input, spatial_coords = spatial_coords)


test_that("nnSVG runs with non-SPE inputs", {
  expect_true(is.matrix(out))
  expect_true(is.numeric(out))
  expect_equal(dim(out), c(10, 14))
})
