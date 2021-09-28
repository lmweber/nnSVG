#' preprocessSVG
#' 
#' Preprocessing steps to prepare input data for nnSVG
#' 
#' Convenience function to run several preprocessing steps to prepare input 
#' data object for nnSVG. This function is designed for data from the 10x 
#' Genomics Visium platform, and is used for examples in the nnSVG package.
#' 
#' In general, the code in this function will need to be adapted for a given 
#' dataset, so we recommend running the steps individually instead of using 
#' this function. The steps are described in detail in our online book 
#' \code{Orchestrating Spatially Resolved Transcriptomics Analysis with 
#' Bioconductor (OSTA)}.
#' 
#' 
#' @param spe \code{SpatialExperiment} Input data, assumed to be a 
#'   \code{\link{SpatialExperiment}} object containing an assay named 
#'   \code{counts} containing raw unique molecular identifier (UMI) counts, and 
#'   spatial coordinates stored in the \code{spatialCoords} slot.
#' 
#' @param in_tissue \code{logical} Whether to keep only spots over tissue, 
#'   identified by the column \code{in_tissue} in the \code{spatialData} slot. 
#'   Default = TRUE.
#' 
#' @param filter_genes \code{integer} Whether to filter low-expressed genes 
#'   according to total unique molecular identifier (UMI) counts per gene 
#'   across spatial coordinates. If an integer value is provided, genes with 
#'   less than or equal to this number of total UMI counts across spatial 
#'   coordinates will be removed. Assumes \code{spe} contains an assay named 
#'   \code{counts} containing raw UMI counts. Default = 20. Set to NULL to 
#'   disable.
#' 
#' @param filter_mito \code{logical} Whether to filter mitochondrial genes. 
#'   Assumes the \code{rowData} slot of \code{spe} contains a column named 
#'   \code{gene_name}, which can be used to identify mitochondrial genes. 
#'   Default = TRUE. Set to FALSE to disable.
#' 
#' 
#' @return Returns a \code{SpatialExperiment} object that can be provided to 
#'   \code{\link{nnSVG}}.
#' 
#' 
#' @importFrom SpatialExperiment spatialData
#' @importFrom SingleCellExperiment counts
#' @importFrom SummarizedExperiment assayNames
#' @importFrom scran quickCluster computeSumFactors
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
#' # set seed for reproducibility
#' set.seed(123)
#' spe <- preprocessSVG(spe)
#' 
#' spe
#' 
preprocessSVG <- function(spe, in_tissue = TRUE, 
                          filter_genes = 20, filter_mito = TRUE) {
  
  stopifnot(isClass(spe, "SpatialExperiment"))
  
  # keep only spots over tissue
  
  if (in_tissue) {
    stopifnot("in_tissue" %in% colnames(spatialData(spe)))
    spe <- spe[, spatialData(spe)$in_tissue == 1]
  }
  
  # gene filtering
  
  if (!is.null(filter_genes) & filter_genes > 0 ) {
    stopifnot("counts" %in% assayNames(spe))
    sums <- rowSums(counts(spe))
    ix_remove <- sums <= filter_genes
    message("removing ", sum(ix_remove), " out of ", nrow(spe), " genes due to low expression")
    spe <- spe[!ix_remove, ]
  }
  
  if (filter_mito) {
    stopifnot("gene_name" %in% colnames(rowData(spe)))
    is_mito <- grepl("(^MT-)|(^mt-)", rowData(spe)$gene_name)
    message("removing ", sum(is_mito), " mitochondrial genes")
    spe <- spe[!is_mito, ]
  }
  
  # normalization and log-transformation
  
  qclus <- quickCluster(spe)
  spe <- computeSumFactors(spe, cluster = qclus)
  spe <- logNormCounts(spe)
  
  # return object
  
  spe
}

