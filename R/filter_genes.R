#' Preprocessing function to filter genes
#' 
#' Preprocessing function to filter low-expressed genes and/or mitochondrial
#' genes for 'nnSVG'.
#' 
#' Preprocessing function to filter low-expressed genes and/or mitochondrial
#' genes for 'nnSVG'.
#' 
#' This function can be used to filter out low-expressed genes and/or
#' mitochondrial genes before additional preprocessing (calculating deviance
#' residuals or log-transformed normalized counts) and running 'nnSVG'.
#' 
#' We use this function in the examples and vignettes in the 'nnSVG' package,
#' and provide default filtering parameter values that are appropriate for 10x
#' Genomics Visium data.
#' 
#' Alternatively, users can also perform filtering and preprocessing separately,
#' and run \code{\link{nnSVG}} on a preprocessed \code{SpatialExperiment}
#' object.
#' 
#' 
#' @param spe \code{SpatialExperiment}: Input data, assumed to be formatted as a
#'   \code{SpatialExperiment} object with an \code{assay} slot named
#'   \code{counts} containing raw expression counts.
#' 
#' @param filter_genes_ncounts \code{numeric}: Filtering parameter for
#'   low-expressed genes. Filtering retains genes containing at least
#'   \code{filter_genes_ncounts} expression counts in at least
#'   \code{filter_genes_pcspots} percent of the total number of spatial
#'   locations (spots). Defaults: \code{filter_genes_ncounts} = 2,
#'   \code{filter_genes_pcspots} = 0.5, i.e. keep genes with at least 2 counts
#'   in at least 0.5 percent of spots. Set to NULL to disable.
#' 
#' @param filter_genes_pcspots \code{numeric}: Second filtering parameter for
#'   low-expressed genes. See \code{filter_genes_ncounts} for details.
#' 
#' @param filter_mito \code{logical}: Whether to filter out mitochondrial genes,
#'   identified by gene names starting with "MT" or "mt". This requires that the
#'   \code{rowData} slot of the input object contains a column named
#'   \code{gene_name}. Default = TRUE. Set to FALSE to disable.
#' 
#' 
#' @return Returns \code{SpatialExperiment} with filtered genes (rows) removed.
#' 
#' 
#' @importFrom SingleCellExperiment counts
#' @importFrom SummarizedExperiment assayNames assays rowData 'rowData<-'
#' 
#' @export
#' 
#' @examples
#' library(SpatialExperiment)
#' library(STexampleData)
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
filter_genes <- function(spe, filter_genes_ncounts = 2, 
                         filter_genes_pcspots = 0.5, 
                         filter_mito = TRUE) {
  
  # filter low-expressed genes
  if (!is.null(filter_genes_ncounts) | !is.null(filter_genes_pcspots)) {
    message("Gene filtering: retaining genes with at least ", 
            filter_genes_ncounts, " counts in at least " , 
            filter_genes_pcspots, " % of spatial locations")
    
    stopifnot("counts" %in% assayNames(spe))
    
    nspots <- ceiling(filter_genes_pcspots / 100 * ncol(spe))
    ix_remove <- rowSums(counts(spe) >= filter_genes_ncounts) < nspots
    message("removed ", sum(ix_remove), " out of ", nrow(spe), 
            " genes due to low expression")
    
    spe <- spe[!ix_remove, ]
  }
  
  # filter mitochondrial genes
  if (filter_mito) {
    message("Gene filtering: removing mitochondrial genes")
    stopifnot("gene_name" %in% colnames(rowData(spe)))
    
    is_mito <- grepl("(^MT-)|(^mt-)", rowData(spe)$gene_name)
    message("removing ", sum(is_mito), " mitochondrial genes")
    
    spe <- spe[!is_mito, ]
  }
  
  spe
}

