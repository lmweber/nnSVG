#' preprocessSVG
#' 
#' Preprocessing steps to prepare input data for nnSVG
#' 
#' Convenience function to run several preprocessing steps to prepare input data
#' object for nnSVG. This function is designed for data from the 10x Genomics
#' Visium platform, and is used for examples in the nnSVG package.
#' 
#' In general, the code in this function will need to be adapted for a given
#' dataset, so we recommend running the steps individually instead of using this
#' function. The steps are described in more detail in our online book
#' "Orchestrating Spatially Resolved Transcriptomics Analysis with Bioconductor
#' (OSTA)".
#' 
#' 
#' @param spe \code{SpatialExperiment} Input data, assumed to be a
#'   \code{\link{SpatialExperiment}} object containing an assay named
#'   \code{counts} containing unique molecular identifier (UMI) counts, and
#'   spatial coordinates stored in the \code{spatialCoords} slot.
#' 
#' @param in_tissue \code{logical} Whether to keep only spots over tissue,
#'   identified by the column \code{in_tissue} in the \code{spatialData} slot.
#'   Default = TRUE.
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
#' @param residuals \code{logical} Whether to calculate deviance residuals from
#'   an approximate multinomial model using the \code{scry} package for input to
#'   BRISC. Default = TRUE.
#' 
#' @param logcounts \code{logical} Whether to calculate log-transformed
#'   normalized counts (logcounts) using the \code{scran} package. Default =
#'   TRUE.
#' 
#' @param deconv \code{logical} Whether to use deconvolution method to calculate
#'   logcounts (see \code{?scran::computeSumFactors}). If \code{FALSE}, library
#'   size normalization will be used instead. Default = TRUE.
#' 
#' 
#' @return Returns a \code{SpatialExperiment} object that can be provided to
#'   \code{\link{nnSVG}}.
#' 
#' 
#' @importFrom SpatialExperiment spatialData 'colData<-'
#' @importFrom SingleCellExperiment counts
#' @importFrom SummarizedExperiment assayNames
#' @importFrom scry nullResiduals
#' @importFrom scran quickCluster computeSumFactors
#' @importFrom scuttle computeLibraryFactors
#' @importFrom scuttle logNormCounts
#' @importFrom methods isClass
#' @importFrom Matrix rowSums
#' 
#' @export
#' 
#' @examples
#' library(SpatialExperiment)
#' library(STexampleData)
#' 
#' spe <- Visium_humanDLPFC()
#' 
#' # subset genes for faster runtime in this example
#' set.seed(123)
#' spe <- spe[sample(seq_len(1000)), ]
#' 
#' # set seed for reproducibility
#' set.seed(123)
#' spe <- preprocessSVG(spe)
#' 
#' spe
#' 
preprocessSVG <- function(spe, in_tissue = TRUE, 
                          filter_genes = 5, filter_mito = TRUE, 
                          residuals = TRUE, logcounts = TRUE, 
                          deconv = TRUE) {
  
  stopifnot(isClass(spe, "SpatialExperiment"))
  
  # keep only spots over tissue
  
  if (in_tissue) {
    stopifnot("in_tissue" %in% colnames(spatialData(spe)))
    spe <- spe[, spatialData(spe)$in_tissue == 1]
  }
  
  # gene filtering
  
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
  
  # calculate deviance residuals or log-transformed normalized counts
  
  if (residuals) {
    spe <- nullResiduals(spe, assay = "counts", 
                         fam = "binomial", type = "deviance")
  }
  if (logcounts) {
    if (deconv) {
      qclus <- quickCluster(spe)
      spe <- computeSumFactors(spe, cluster = qclus)
      spe <- logNormCounts(spe)
    } else {
      spe <- computeLibraryFactors(spe)
      spe <- logNormCounts(spe)
    }
  }
  
  # return object
  
  spe
}

