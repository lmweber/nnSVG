#' nnSVG
#' 
#' Run 'nnSVG' to identify spatially variable genes
#' 
#' Run 'nnSVG' method to identify spatially variable genes (SVGs) in spatially
#' resolved transcriptomics data.
#' 
#' The method is based on nearest neighbor Gaussian processes (Datta et al.
#' 2016) and uses the BRISC algorithm (Saha and Datta 2018). The method scales
#' linearly in the number of spatial coordinates, and can be applied to datasets
#' containing thousands of spatial coordinates.
#' 
#' This function runs BRISC once per gene for model fitting and parameter
#' estimation, using parallelization for faster runtime with one core per BRISC
#' run. The spatial covariance parameter estimates (sigma.sq, tau.sq, phi) from
#' BRISC are stored in 'Theta' in the BRISC output.
#' 
#' nnSVG performs inference on the 'sigma.sq' estimates using an approximate
#' likelihood ratio test against a model without spatial terms, and uses these
#' likelihood ratios to rank SVGs. We also calculate an effect size, defined as
#' the proportion of spatial variance out of total variance, i.e. 'prop_sv =
#' sigma.sq / (sigma.sq + tau.sq)'.
#' 
#' Likelihood ratio tests are calculated using the asymptotic chi-squared
#' distribution with 2 degrees of freedom, and adjusted p-values using the
#' Benjamini-Hochberg method.
#' 
#' Assumes the input is provided as a \code{\link{SpatialExperiment}} object
#' with an \code{assay} slot containing either deviance residuals (e.g. from the
#' \code{scry} package) or log-transformed normalized counts (e.g. from the
#' \code{scran} package), and which has been filtered to remove low-quality
#' spatial coordinates.
#' 
#' Low-expressed genes can be filtered before providing the input to
#' \code{nnSVG()}, or using the default filtering arguments in \code{nnSVG()}.
#' 
#' 
#' @param spe \code{SpatialExperiment} Input data, assumed to be a
#'   \code{\link{SpatialExperiment}} object with \code{assay} slots containing
#'   deviance residuals and/or log-transformed normalized counts, and spatial
#'   coordinates stored in the \code{spatialCoords} slot.
#' 
#' @param x \code{numeric matrix} Optional matrix of covariates (e.g. known cell
#'   types) per spatial coordinate. Number of rows must match the number of
#'   spatial coordinates (columns) in the input object \code{spe}. Default =
#'   NULL, which fits an intercept-only model.
#' 
#' @param assay_name \code{character} Name of assay containing preprocessed
#'   expression values to use for clustering, i.e. either deviance residuals or
#'   log-transformed normalized counts. Assumed to be either
#'   \code{binomial_deviance_residuals} or \code{logcounts}. Default =
#'   \code{binomial_deviance_residuals}.
#' 
#' @param filter_genes \code{integer} Whether to filter low-expressed genes. If
#'   a value is provided, genes with at least 1 unique molecular identifier
#'   (UMI) count in at least this percentage of spatial coordinates will be
#'   kept. Assumes \code{spe} contains an assay named \code{counts} containing
#'   UMI counts. Default = 5, i.e. keep genes with at least 1 UMI in 5% of
#'   spatial coordinates. Set to NULL to disable.
#' 
#' @param filter_mito \code{logical} Whether to filter mitochondrial genes.
#'   Assumes the \code{rowData} slot of \code{spe} contains a column named
#'   \code{gene_name}, which can be used to identify mitochondrial genes.
#'   Default = TRUE. Set to FALSE to disable.
#' 
#' @param n_threads \code{integer} Number of threads for parallelization.
#'   Default = 1.
#' 
#' @param verbose \code{logical} Whether to display verbose output from
#'   \code{BRISC}. Default = FALSE.
#' 
#' 
#' @return Returns output values (including parameter estimates, likelihood
#'   ratios, and adjusted p-values) in \code{rowData} in the \code{spe}
#'   \code{SpatialExperiment} object.
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
#' # filter low-expressed genes
#' filter_genes <- 5
#' n_spots <- ceiling(filter_genes / 100 * ncol(spe))
#' ix_remove <- rowSums(counts(spe) > 0) < n_spots
#' spe <- spe[!ix_remove, ]
#' 
#' # filter mitochondrial genes
#' is_mito <- grepl("(^MT-)|(^mt-)", rowData(spe)$gene_name)
#' spe <- spe[!is_mito, ]
#' 
#' # calculate deviance residuals using scry package
#' set.seed(123)
#' spe <- nullResiduals(spe, assay = "counts", 
#'                      fam = "binomial", type = "deviance")
#' 
#' # subset small number of genes for faster runtime in this example
#' set.seed(123)
#' spe <- spe[sample(seq_len(6)), ]
#' 
#' # run nnSVG
#' # note: gene filtering was already performed above, so disable it here
#' set.seed(123)
#' spe <- nnSVG(spe, filter_genes = FALSE, filter_mito = FALSE, n_threads = 1)
#' 
#' # show results
#' # for more details see extended example in vignette
#' rowData(spe)
#' 
nnSVG <- function(spe, x = NULL, 
                  assay_name = c("binomial_deviance_residuals", "logcounts"), 
                  filter_genes = 5, filter_mito = TRUE, 
                  n_threads = 1, verbose = FALSE) {
  
  if (!is.null(x)) stopifnot(nrow(x) == ncol(spe))
  
  assay_name <- match.arg(assay_name, c("binomial_deviance_residuals", "logcounts"))
  stopifnot(assay_name %in% assayNames(spe))
  
  # --------------
  # gene filtering
  # --------------
  
  if (!is.null(filter_genes) & filter_genes > 0) {
    stopifnot("counts" %in% assayNames(spe))
    n_spots <- ceiling(filter_genes / 100 * ncol(spe))
    ix_remove <- rowSums(counts(spe) > 0) < n_spots
    message("removing ", sum(ix_remove), " out of ", nrow(spe), " genes due to low expression")
    spe <- spe[!ix_remove, ]
  }
  
  if (filter_mito) {
    stopifnot("gene_name" %in% colnames(rowData(spe)))
    is_mito <- grepl("(^MT-)|(^mt-)", rowData(spe)$gene_name)
    message("removing ", sum(is_mito), " mitochondrial genes")
    spe <- spe[!is_mito, ]
  }
  
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
        out_i <- BRISC_estimation(coords = coords, y = y_i, x = x, 
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
      mat_brisc, 
      mean = NA, 
      var = NA, 
      spcov = NA
    )
  }
  
  # ratio of spatial to non-spatial variance
  mat_brisc <- cbind(
    mat_brisc, 
    ratio_sv = mat_brisc[, "sigma.sq"] / mat_brisc[, "tau.sq"]
  )
  
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
    if (is.null(x)) {
      x <- rep(1, ncol(spe))
    }
    as.numeric(logLik(lm(y_i ~ x)))
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

