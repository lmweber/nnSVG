#' nnSVG
#' 
#' Function to run 'nnSVG' method to identify spatially variable genes (SVGs) in
#' spatially-resolved transcriptomics data.
#' 
#' Function to run 'nnSVG' method to identify spatially variable genes (SVGs) in
#' spatially-resolved transcriptomics data.
#' 
#' The 'nnSVG' method is based on nearest-neighbor Gaussian processes (Datta et
#' al. 2016) and uses the BRISC algorithm (Saha and Datta 2018) for model
#' fitting and parameter estimation. The method scales linearly with the number
#' of spatial locations, and can be applied to datasets containing thousands or
#' more spatial locations. For more details on the method, see our paper.
#' 
#' This function runs 'nnSVG' for a full dataset. The function fits a separate
#' model for each gene, using parallelization with BiocParallel for faster
#' runtime. The parameter estimates from BRISC (sigma.sq, tau.sq, phi) for each
#' gene are stored in 'Theta' in the BRISC output.
#' 
#' Note that the method and this function are designed for a single tissue
#' section. For an example of how to run nnSVG in a dataset consisting of
#' multiple tissue sections, see the tutorial in the nnSVG package vignette.
#' 
#' 'nnSVG' performs inference on the spatial variance parameter estimates
#' (sigma.sq) using a likelihood ratio (LR) test against a simpler linear model
#' without spatial terms (i.e. without tau.sq or phi). The estimated LR
#' statistics can then be used to rank SVGs. P-values are calculated from the LR
#' statistics using the asymptotic chi-squared distribution with 2 degrees of
#' freedom, and multiple testing adjusted p-values are calculated using the
#' Benjamini-Hochberg method. We also calculate an effect size, defined as the
#' proportion of spatial variance, 'prop_sv = sigma.sq / (sigma.sq + tau.sq)'.
#' 
#' The function assumes the input is provided either as a
#' \code{SpatialExperiment} object or a \code{numeric} matrix of values. If the
#' input is a \code{SpatialExperiment} object, it is assumed to contain an
#' \code{assay} slot containing either log-transformed normalized counts (also
#' known as logcounts, e.g. from the \code{scran} package) or deviance residuals
#' (e.g. from the \code{scry} package), which have been preprocessed, quality
#' controlled, and filtered to remove low-quality spatial locations. If the
#' input is a \code{numeric} matrix of values, these values are assumed to
#' already be normalized and transformed (e.g. logcounts).
#' 
#' 
#' @param input \code{SpatialExperiment} or \code{numeric} matrix: Input data,
#'   which can either be a \code{SpatialExperiment} object or a \code{numeric}
#'   matrix of values. If it is a \code{SpatialExperiment} object, it is assumed
#'   to have an \code{assay} slot containing either logcounts (e.g. from the
#'   \code{scran} package) or deviance residuals (e.g. from the \code{scry}
#'   package), and a \code{spatialCoords} slot containing spatial coordinates of
#'   the measurements. If it is a \code{numeric} matrix, the values are assumed
#'   to already be normalized and transformed (e.g. logcounts), formatted as
#'   \code{rows = genes} and \code{columns = spots}, and a separate
#'   \code{numeric} matrix of spatial coordinates must also be provided with the
#'   \code{spatial_coords} argument.
#' 
#' @param spatial_coords \code{numeric} matrix: Matrix containing columns of
#'   spatial coordinates, formatted as \code{rows = spots}. This must be
#'   provided if \code{input} is provied as a \code{numeric} matrix of values,
#'   and is ignored if \code{input} is provided as a \code{SpatialExperiment}
#'   object. Default = NULL.
#' 
#' @param X \code{numeric} matrix: Optional design matrix containing columns of
#'   covariates per spatial location, e.g. known spatial domains. Number of rows
#'   must match the number of spatial locations. Default = NULL, which fits an
#'   intercept-only model.
#' 
#' @param assay_name \code{character}: If \code{input} is provided as a
#'   \code{SpatialExperiment} object, this argument selects the name of the
#'   \code{assay} slot in the input object containing the preprocessed gene
#'   expression values. For example, \code{logcounts} for log-transformed
#'   normalized counts from the \code{scran} package, or
#'   \code{binomial_deviance_residuals} for deviance residuals from the
#'   \code{scry} package. Default = \code{"logcounts"}, or ignored if
#'   \code{input} is provided as a \code{numeric} matrix of values.
#' 
#' @param n_neighbors \code{integer}: Number of nearest neighbors for fitting
#'   the nearest-neighbor Gaussian process (NNGP) model with BRISC. The default
#'   value is 10, which we recommend for most datasets. Higher numbers (e.g. 15)
#'   may give slightly improved likelihood estimates in some datasets (at the
#'   expense of slower runtime), and smaller numbers (e.g. 5) will give faster
#'   runtime (at the expense of reduced performance). Default = 10.
#' 
#' @param order \code{character}: Ordering scheme to use for ordering
#'   coordinates with BRISC. Default = \code{"AMMD"} for "approximate maximum
#'   minimum distance", which is recommended for datasets with at least 65
#'   spots. For very small datasets (n <= 65), \code{"Sum_coords"} can be used
#'   instead. See BRISC documentation for details. Default = \code{"AMMD"}.
#' 
#' @param n_threads \code{integer}: Number of threads for parallelization.
#'   Default = 1. We recommend setting this equal to the number of cores
#'   available (if working on a laptop or desktop) or around 10 or more (if
#'   working on a compute cluster).
#' 
#' @param BPPARAM \code{BiocParallelParam}: Optional additional argument for
#'   parallelization. This argument is provided for advanced users of
#'   \code{BiocParallel} for further flexibility for parallelization on some
#'   operating systems. If provided, this should be an instance of
#'   \code{BiocParallelParam}. For most users, the recommended option is to use
#'   the \code{n_threads} argument instead. Default = NULL, in which case
#'   \code{n_threads} will be used instead.
#' 
#' @param verbose \code{logical}: Whether to display verbose output for model
#'   fitting and parameter estimation from \code{BRISC}. Default = FALSE.
#' 
#' 
#' @return If the input was provided as a \code{SpatialExperiment} object, the
#'   output values are returned as additional columns in the \code{rowData} slot
#'   of the input object. If the input was provided as a \code{numeric} matrix
#'   of values, the output is returned as a \code{numeric} matrix. The output
#'   values include spatial variance parameter estimates, likelihood ratio (LR)
#'   statistics, effect sizes (proportion of spatial variance), p-values, and
#'   multiple testing adjusted p-values.
#' 
#' 
#' @importFrom SpatialExperiment spatialCoords
#' @importFrom SingleCellExperiment counts logcounts
#' @importFrom SummarizedExperiment assayNames assays rowData 'rowData<-'
#' @importFrom BRISC BRISC_order BRISC_neighbor BRISC_estimation
#' @importFrom BiocParallel bplapply MulticoreParam
#' @importFrom Matrix rowSums rowMeans colSums
#' @importFrom matrixStats rowVars
#' @importFrom stats lm logLik pchisq p.adjust
#' @importFrom methods is
#' 
#' @export
#' 
#' @examples
#' library(SpatialExperiment)
#' library(STexampleData)
#' library(scran)
#' 
#' 
#' ### Example 1
#' ### for more details see extended example in vignette
#' 
#' # load example dataset from STexampleData package
#' spe1 <- Visium_humanDLPFC()
#' 
#' # preprocessing steps
#' # keep only spots over tissue
#' spe1 <- spe1[, colData(spe1)$in_tissue == 1]
#' # skip spot-level quality control (already performed in this dataset)
#' # filter low-expressed and mitochondrial genes
#' spe1 <- filter_genes(spe1)
#' # calculate logcounts using library size factors
#' spe1 <- computeLibraryFactors(spe1)
#' spe1 <- logNormCounts(spe1)
#' 
#' # select small number of genes for fast runtime in this example
#' set.seed(123)
#' ix <- c(
#'   which(rowData(spe1)$gene_name %in% c("PCP4", "NPY")), 
#'   sample(seq_len(nrow(spe1)), 2)
#' )
#' spe1 <- spe1[ix, ]
#' 
#' # run nnSVG
#' set.seed(123)
#' spe1 <- nnSVG(spe1)
#' 
#' # show results
#' rowData(spe1)
#' 
#' 
#' ### Example 2: With covariates
#' ### for more details see extended example in vignette
#' 
#' # load example dataset from STexampleData package
#' spe2 <- SlideSeqV2_mouseHPC()
#' 
#' # preprocessing steps
#' # remove spots with NA cell type labels
#' spe2 <- spe2[, !is.na(colData(spe2)$celltype)]
#' # skip spot-level quality control (already performed in this dataset)
#' # filter low-expressed and mitochondrial genes
#' spe2 <- filter_genes(
#'   spe2, filter_genes_ncounts = 1, filter_genes_pcspots = 1, 
#'   filter_mito = TRUE
#' )
#' # calculate logcounts using library size normalization
#' spe2 <- computeLibraryFactors(spe2)
#' spe2 <- logNormCounts(spe2)
#' 
#' # select small number of genes for fast runtime in this example
#' set.seed(123)
#' ix <- c(
#'   which(rowData(spe2)$gene_name %in% c("Cpne9", "Rgs14")), 
#'   sample(seq_len(nrow(spe2)), 2)
#' )
#' spe2 <- spe2[ix, ]
#' 
#' # create model matrix for cell type labels
#' X <- model.matrix(~ colData(spe2)$celltype)
#' 
#' # run nnSVG with covariates
#' set.seed(123)
#' spe2 <- nnSVG(spe2, X = X)
#' 
#' # show results
#' rowData(spe2)
#' 
nnSVG <- function(input, spatial_coords = NULL, X = NULL, 
                  assay_name = "logcounts", 
                  n_neighbors = 10, order = "AMMD", 
                  n_threads = 1, BPPARAM = NULL, 
                  verbose = FALSE) {
  
  if (is(input, "SpatialExperiment")) {
    spe <- input
    stopifnot(assay_name %in% assayNames(spe))
    # check for rows or columns of all zero counts
    if ("counts" %in% assayNames(spe)) {
      if (sum(rowSums(counts(spe)) == 0) > 0 | sum(colSums(counts(spe)) == 0) > 0) {
        warning("Rows (genes) and/or columns (spots) containing all zero counts ", 
                "have been found. Please see examples in tutorial for code to ", 
                "filter out zeros and/or low-expressed genes to avoid errors.")
      }
    }
  }
  
  if (!is.null(X)) {
    stopifnot(nrow(X) == ncol(input))
  }
  
  if (is.null(BPPARAM)) {
    BPPARAM <- MulticoreParam(workers = n_threads)
  }
  
  # -----------------------
  # run BRISC for each gene
  # -----------------------
  
  if (is(input, "SpatialExperiment")) {
    y <- assays(spe)[[assay_name]]
    coords <- spatialCoords(spe)
  } else {
    y <- input
    coords <- spatial_coords
    row_names <- rownames(input)
  }
  
  # scale coordinates proportionally
  range_all <- max(apply(coords, 2, function(col) diff(range(col))))
  coords <- apply(coords, 2, function(col) (col - min(col)) / range_all)
  
  # calculate ordering of coordinates
  order_brisc <- BRISC_order(coords, order = order, verbose = verbose)
  
  # calculate nearest neighbors
  nn_brisc <- BRISC_neighbor(coords, n.neighbors = n_neighbors, n_omp = 1, 
                             search.type = "tree", ordering = order_brisc, 
                             verbose = verbose)
  
  # run BRISC using parallelization
  ix <- seq_len(nrow(y))
  out_brisc <- bplapply(ix, function(i) {
    # fit model (intercept-only model if x is NULL)
    y_i <- y[i, ]
    suppressWarnings({
      runtime <- system.time({
        out_i <- BRISC_estimation(coords = coords, y = y_i, x = X, 
                                  cov.model = "exponential", 
                                  ordering = order_brisc, neighbor = nn_brisc, 
                                  verbose = verbose)
      })
    })
    res_i <- c(
      out_i$Theta, 
      loglik = out_i$log_likelihood, 
      runtime = runtime[["elapsed"]]
    )
    res_i
  }, BPPARAM = BPPARAM)
  
  # collapse output list into matrix
  mat_brisc <- do.call("rbind", out_brisc)
  
  # --------------------
  # calculate statistics
  # --------------------
  
  if (is(input, "SpatialExperiment") && ("logcounts" %in% assayNames(spe))) {
    lc <- logcounts(spe)
    # mean logcounts
    mat_brisc <- cbind(
      mat_brisc, 
      mean = rowMeans(lc)
    )
    # variance of logcounts
    mat_brisc <- cbind(
      mat_brisc, 
      var = rowVars(as.matrix(lc))
    )
    # spatial coefficient of variation
    mat_brisc <- cbind(
      mat_brisc, 
      spcov = sqrt(mat_brisc[, "sigma.sq"]) / mat_brisc[, "mean"]
    )
  } else {
    # return NAs if logcounts not provided
    mat_brisc <- cbind(
      mat_brisc, mean = NA, var = NA, spcov = NA
    )
  }
  
  # proportion of spatial variance out of total variance
  mat_brisc <- cbind(
    mat_brisc, 
    prop_sv = mat_brisc[, "sigma.sq"] / (mat_brisc[, "sigma.sq"] + mat_brisc[, "tau.sq"])
  )
  
  # ------------------------------------------
  # likelihood ratio (LR) statistics and tests
  # ------------------------------------------
  
  if (is(input, "SpatialExperiment")) {
    nrows <- nrow(spe)
    ncols <- ncol(spe)
  } else {
    nrows <- nrow(input)
    ncols <- ncol(input)
  }
  
  # calculate log likelihoods for nonspatial models
  
  loglik_lm <- vapply(seq_len(nrows), function(i) {
    y_i <- y[i, ]
    if (is.null(X)) {
      X <- rep(1, ncols)
    }
    # model formula without intercept to enable weighted model
    as.numeric(logLik(lm(y_i ~ X - 1)))
  }, numeric(1))
  
  mat_brisc <- cbind(
    mat_brisc, 
    loglik_lm = loglik_lm
  )
  
  # calculate LR statistics and tests (Wilks' theorem, asymptotic chi-square
  # with 2 degrees of freedom)
  
  LR_stat <- -2 * (mat_brisc[, "loglik_lm"] - mat_brisc[, "loglik"])
  
  pval <- 1 - pchisq(LR_stat, df = 2)
  padj <- p.adjust(pval, method = "BH")
  
  # rank SVGs according to LR statistics
  LR_rank <- rank(-1 * LR_stat)
  
  mat_brisc <- cbind(
    mat_brisc, 
    LR_stat = LR_stat, 
    rank = LR_rank, 
    pval = pval, 
    padj = padj
  )
  
  # --------------
  # return outputs
  # --------------
  
  if (is(input, "SpatialExperiment")) {
    # return in rowData of spe object
    stopifnot(nrow(spe) == nrow(mat_brisc))
    rowData(spe) <- cbind(rowData(spe), mat_brisc)
    spe
  } else {
    # return as numeric matrix
    stopifnot(nrow(input) == nrow(mat_brisc))
    rownames(mat_brisc) <- row_names
    mat_brisc
  }
}

