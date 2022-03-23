#' nnSVG
#' 
#' Function to run 'nnSVG' method to identify spatially variable genes (SVGs) in
#' spatially resolved transcriptomics (ST) data.
#' 
#' Function to run 'nnSVG' method to identify spatially variable genes (SVGs) in
#' spatially resolved transcriptomics (ST) data.
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
#' 'nnSVG' performs inference on the spatial variance parameter estimates
#' (sigma.sq) using a likelihood ratio (LR) test against a simpler linear model
#' without spatial terms (i.e. without tau.sq or phi). The estimated LR
#' statistics can then be used to rank SVGs. P-values are calculated from the LR
#' statistics using the asymptotic chi-squared distribution with 2 degrees of
#' freedom, and multiple testing adjusted p-values are calculated using the
#' Benjamini-Hochberg method. We also calculate an effect size, defined as the
#' proportion of spatial variance, 'prop_sv = sigma.sq / (sigma.sq + tau.sq)'.
#' 
#' The function assumes the input is provided as a \code{SpatialExperiment}
#' object containing an \code{assay} slot containing either log-transformed
#' normalized counts (e.g. from the \code{scran} package) or deviance residuals
#' (e.g. from the \code{scry} package), which have been preprocessed, quality
#' controlled, and filtered to remove low-quality spatial locations.
#' 
#' 
#' @param spe \code{SpatialExperiment}: Input data, assumed to be formatted as a
#'   \code{SpatialExperiment} object with an \code{assay} slot containing either
#'   log-transformed normalized counts (e.g. from the \code{scran} package) or
#'   deviance residuals (e.g. from the \code{scry} package), and a
#'   \code{spatialCoords} slot containing spatial coordinates of the
#'   measurements.
#' 
#' @param X \code{numeric} matrix: Optional design matrix containing columns of
#'   covariates per spatial location, e.g. known spatial domains. Number of rows
#'   must match the number of spatial locations. Default = NULL, which fits an
#'   intercept-only model.
#' 
#' @param assay_name \code{character}: Name of the \code{assay} slot in the
#'   input object containing the preprocessed gene expression values. For
#'   example, \code{logcounts} for log-transformed normalized counts from the
#'   \code{scran} package, or \code{binomial_deviance_residuals} for deviance
#'   residuals from the \code{scry} package. Default = \code{logcounts}.
#' 
#' @param n_neighbors \code{integer}: Number of nearest neighbors for fitting
#'   the nearest-neighbor Gaussian process (NNGP) model with BRISC. The default
#'   value is 15, which works well in most settings. Smaller numbers (e.g. 5-10)
#'   will give faster runtime at the expense of reduced performance. Default =
#'   15.
#' 
#' @param n_threads \code{integer}: Number of threads for parallelization.
#'   Default = 1.
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
#' @return Returns output values as additional columns in the \code{rowData}
#'   slot of the input object, including spatial variance parameter estimates,
#'   likelihood ratio (LR) statistics, effect sizes (proportion of spatial
#'   variance), p-values, and multiple testing adjusted p-values.
#' 
#' 
#' @importFrom SpatialExperiment spatialCoords
#' @importFrom SingleCellExperiment counts logcounts
#' @importFrom SummarizedExperiment assayNames assays rowData 'rowData<-'
#' @importFrom BRISC BRISC_order BRISC_neighbor BRISC_estimation
#' @importFrom BiocParallel bplapply MulticoreParam
#' @importFrom Matrix rowSums rowMeans
#' @importFrom matrixStats rowVars
#' @importFrom stats lm logLik pchisq p.adjust
#' 
#' @export
#' 
#' @examples
#' library(SpatialExperiment)
#' library(STexampleData)
#' library(scran)
#' 
#' # load example dataset from STexampleData package
#' spe <- Visium_humanDLPFC()
#' 
#' # preprocessing steps
#' 
#' # keep only spots over tissue
#' spe <- spe[, colData(spe)$in_tissue == 1]
#' 
#' # skip spot-level quality control, since this has been performed previously 
#' # on this dataset
#' 
#' # filter low-expressed and mitochondrial genes
#' spe <- filter_genes(spe)
#' 
#' # calculate log-transformed normalized counts using scran package
#' set.seed(123)
#' qclus <- quickCluster(spe)
#' spe <- computeSumFactors(spe, cluster = qclus)
#' spe <- logNormCounts(spe)
#' 
#' # select small number of genes for faster runtime in this example
#' set.seed(123)
#' ix <- sample(seq_len(nrow(spe)), 4)
#' spe <- spe[ix, ]
#' 
#' # run nnSVG
#' set.seed(123)
#' spe <- nnSVG(spe)
#' 
#' # show results
#' # for more details see extended example in vignette
#' rowData(spe)
#' 
nnSVG <- function(spe, X = NULL, 
                  assay_name = "logcounts", n_neighbors = 15, 
                  n_threads = 1, BPPARAM = NULL, 
                  verbose = FALSE) {
  
  if (!is.null(X)) stopifnot(nrow(X) == ncol(spe))
  
  stopifnot(assay_name %in% assayNames(spe))
  
  if (is.null(BPPARAM)) {
    BPPARAM <- MulticoreParam(workers = n_threads)
  }
  
  # -----------------------
  # run BRISC for each gene
  # -----------------------
  
  y <- assays(spe)[[assay_name]]
  
  # scale coordinates proportionally
  coords <- spatialCoords(spe)
  range_all <- max(apply(coords, 2, function(col) diff(range(col))))
  coords <- apply(coords, 2, function(col) (col - min(col)) / range_all)
  
  # calculate ordering of coordinates
  order_brisc <- BRISC_order(coords, order = "Sum_coords", verbose = verbose)
  
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
  
  if ("logcounts" %in% assayNames(spe)) {
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
  
  # calculate log likelihoods for nonspatial models
  
  loglik_lm <- vapply(seq_len(nrow(spe)), function(i) {
    y_i <- y[i, ]
    if (is.null(X)) {
      X <- rep(1, ncol(spe))
    }
    as.numeric(logLik(lm(y_i ~ X)))
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
  
  # -------------------------------
  # return in rowData of spe object
  # -------------------------------
  
  stopifnot(nrow(spe) == nrow(mat_brisc))
  
  rowData(spe) <- cbind(rowData(spe), mat_brisc)
  
  spe
}

