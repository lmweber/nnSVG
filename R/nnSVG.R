#' nnSVG
#' 
#' Function to run 'nnSVG' method to identify spatially variable genes (SVGs) in
#' spatially resolved transcriptomics (ST) data.
#' 
#' Function to run 'nnSVG' method to identify spatially variable genes (SVGs) in
#' spatially resolved transcriptomics (ST) data.
#' 
#' The 'nnSVG' method is based on nearest-neighbor Gaussian processes (Datta et
#' al. 2016) and uses the BRISC algorithm (Saha and Datta 2018). The method
#' scales linearly in the number of spatial locations, and can be applied to
#' datasets containing thousands or more spatial locations. For more details on
#' the method, see our paper.
#' 
#' This function runs 'nnSVG' for a full dataset. Internally, the function calls
#' BRISC once per gene to perform model fitting and parameter estimation, using
#' parallelization with one core per BRISC run. The spatial covariance parameter
#' estimates (sigma.sq, tau.sq, phi) are stored in 'Theta' in the BRISC output.
#' 
#' 'nnSVG' then performs inference on the 'sigma.sq' estimates using an
#' approximate likelihood ratio (LR) test against a model without spatial terms,
#' and uses the estimated LR statistics to rank SVGs. We also calculate an
#' effect size, defined as the proportion of spatial variance, 'prop_sv =
#' sigma.sq / (sigma.sq + tau.sq)'.
#' 
#' LR tests are calculated using the asymptotic chi-squared distribution with 2
#' degrees of freedom. Multiple testing adjusted p-values are calculated using
#' the Benjamini-Hochberg method.
#' 
#' The function assumes the input is provided as a \code{SpatialExperiment}
#' object containing an \code{assay} slot containing either deviance residuals
#' (e.g. from the \code{scry} package) or log-transformed normalized counts
#' (e.g. from the \code{scran} package), which have been preprocessed, quality
#' controlled, and filtered to remove any low-quality spatial locations.
#' 
#' 
#' @param spe \code{SpatialExperiment}: Input data, assumed to be formatted as a
#'   \code{SpatialExperiment} object with an \code{assay} slot containing either
#'   deviance residuals (e.g. from the \code{scry} package) or log-transformed
#'   normalized counts (e.g. from the \code{scran} package), and a
#'   \code{spatialCoords} slot containing spatial coordinates of the
#'   measurements.
#' 
#' @param X \code{numeric matrix}: Optional design matrix containing columns of
#'   covariates per spatial location, e.g. known spatial domains. Number of rows
#'   must match the number of spatial locations. Default = NULL, which fits an
#'   intercept-only model.
#' 
#' @param assay_name \code{character}: Name of the \code{assay} slot in the
#'   input object containing the preprocessed gene expression values. For
#'   example, \code{binomial_deviance_residuals} for deviance residuals from the
#'   \code{scry} package, or \code{logcounts} for log-transformed normalized
#'   counts from the \code{scran} package. Default =
#'   \code{binomial_deviance_residuals}.
#' 
#' @param n_threads \code{integer}: Number of threads for parallelization.
#'   Default = 1.
#' 
#' @param verbose \code{logical}: Whether to display verbose output for model
#'   fitting and parameter estimation from \code{BRISC}. Default = FALSE.
#' 
#' 
#' @return Returns output values as additional columns in the \code{rowData}
#'   slot of the input object, including spatial covariance parameter estimates,
#'   likelihood ratio (LR) statistics, effect sizes (proportion of spatial
#'   covariance), p-values, and multiple testing adjusted p-values.
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
#' library(scry)
#' 
#' # load example dataset from STexampleData package
#' spe <- Visium_humanDLPFC()
#' 
#' # preprocessing steps
#' 
#' # keep only spots over tissue
#' spe <- spe[, colData(spe)$in_tissue == 1]
#' 
#' # filter low-expressed and mitochondrial genes
#' spe <- filter_genes(spe)
#' 
#' # calculate deviance residuals using scry package
#' spe <- nullResiduals(spe, assay = "counts", 
#'                      fam = "binomial", type = "deviance")
#' 
#' # select small set of random genes and several known SVGs for faster runtime 
#' # in this example
#' set.seed(123)
#' ix_20 <- sample(seq_len(nrow(spe)), 20)
#' known_genes <- c("MOBP", "PCP4", "SNAP25", "HBB", "IGKC", "NPY")
#' ix_known <- which(rowData(spe)$gene_name %in% known_genes)
#' ix <- c(ix_20, ix_known)
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
                  assay_name = "binomial_deviance_residuals", 
                  n_threads = 1, verbose = FALSE) {
  
  if (!is.null(X)) stopifnot(nrow(X) == ncol(spe))
  
  stopifnot(assay_name %in% assayNames(spe))
  
  # ---------
  # run BRISC
  # ---------
  
  y <- assays(spe)[[assay_name]]
  
  # scale coordinates proportionally
  coords <- spatialCoords(spe)
  range_all <- max(apply(coords, 2, function(col) diff(range(col))))
  coords <- apply(coords, 2, function(col) (col - min(col)) / range_all)
  
  # calculate ordering of coordinates
  order_brisc <- BRISC_order(coords, order = "AMMD", verbose = verbose)
  
  # calculate nearest neighbors
  nn_brisc <- BRISC_neighbor(coords, n.neighbors = 15, n_omp = 1, 
                             search.type = "cb", ordering = order_brisc, 
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
  }, BPPARAM = MulticoreParam(workers = n_threads))
  
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
  
  # calculate log likelihoods for non-spatial models
  
  loglik_lm <- sapply(seq_len(nrow(spe)), function(i) {
    y_i <- y[i, ]
    if (is.null(X)) {
      X <- rep(1, ncol(spe))
    }
    as.numeric(logLik(lm(y_i ~ X)))
  })
  
  mat_brisc <- cbind(
    mat_brisc, 
    loglik_lm = loglik_lm
  )
  
  # calculate LR statistics and tests (Wilks' theorem, asymptotic chi-square 
  # with 2 degrees of freedom since 2 more parameters in full model)
  
  LR_stat = -2 * (mat_brisc[, "loglik_lm"] - mat_brisc[, "loglik"])
  
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

