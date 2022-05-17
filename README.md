# nnSVG

[![R build status](https://github.com/lmweber/nnSVG/workflows/R-CMD-check-bioc/badge.svg)](https://github.com/lmweber/nnSVG/actions)


## Overview

`nnSVG` is a method for scalable identification of spatially variable genes (SVGs) in spatially-resolved transcriptomics data.

The `nnSVG` method is based on nearest-neighbor Gaussian processes ([Datta et al., 2016](https://www.tandfonline.com/doi/full/10.1080/01621459.2015.1044091), [Finley et al., 2019](https://www.tandfonline.com/doi/full/10.1080/10618600.2018.1537924)) and uses the BRISC algorithm ([Saha and Datta, 2018](https://onlinelibrary.wiley.com/doi/full/10.1002/sta4.184)) for model fitting and parameter estimation. `nnSVG` allows identification and ranking of SVGs with flexible length scales across a tissue slide or within spatial domains defined by covariates. The method scales linearly with the number of spatial locations and can be applied to datasets containing thousands or more spatial locations.

`nnSVG` is implemented as an R package within the Bioconductor framework, and is available from [Bioconductor](https://bioconductor.org/packages/nnSVG).

Our paper describing the method is available from [bioRxiv](https://www.biorxiv.org/content/10.1101/2022.05.16.492124v1).


## Installation

The package can be installed from [Bioconductor](https://bioconductor.org/packages/nnSVG) as follows, using R version 4.2 or above:

```
install.packages("BiocManager")
BiocManager::install("nnSVG")
```

The latest development version of the package can also be installed from GitHub:

```
remotes::install_github("lmweber/nnSVG")
```


## Workflow

An example workflow is available in the package vignette from [Bioconductor](https://bioconductor.org/packages/nnSVG). A short overview is also provided below.

### Load packages

```
library(nnSVG)
library(STexampleData)
library(scran)
library(ggplot2)
```

### Load example dataset

```
# load example dataset from STexampleData package
spe <- Visium_humanDLPFC()
dim(spe)
```

### Preprocessing

```
# keep spots over tissue
spe <- spe[, colData(spe)$in_tissue == 1]
dim(spe)

# spot-level quality control: already performed on this example dataset

# filter low-expressed and mitochondrial genes
# using function from nnSVG package with default filtering parameters
spe <- filter_genes(spe)
dim(spe)

# calculate log-transformed normalized counts (logcounts) using scran package
set.seed(123)
qclus <- quickCluster(spe)
spe <- computeSumFactors(spe, cluster = qclus)
spe <- logNormCounts(spe)
assayNames(spe)
```

### Subset dataset for this example

```
# select small set of random genes and several known SVGs for faster runtime in this example workflow
set.seed(123)
ix_random <- sample(seq_len(nrow(spe)), 10)
known_genes <- c("MOBP", "PCP4", "SNAP25", "HBB", "IGKC", "NPY")
ix_known <- which(rowData(spe)$gene_name %in% known_genes)
ix <- c(ix_known, ix_random)

spe <- spe[ix, ]
dim(spe)
```

### Run nnSVG

```
# set seed for reproducibility
set.seed(123)
# run nnSVG using a single thread for this example workflow
spe <- nnSVG(spe, n_threads = 1)

# show results
rowData(spe)
```

### Investigate results

The results are stored in the `rowData` of the `SpatialExperiment` object.

The main results of interest are:

- `LR_stat`: likelihood ratio (LR) statistics used to rank SVGs
- `rank`: rank of top SVGs according to LR statistics
- `pval`: approximate p-values
- `padj`: approximate p-values adjusted for multiple testing
- `prop_sv`: effect size defined as proportion of spatial variance

```
# number of significant SVGs
table(rowData(spe)$padj <= 0.05)

# show results for top n SVGs
n <- 10
rowData(spe)[order(rowData(spe)$rank)[1:n], ]
```

### Plot expression

Plot expression of the top-ranked SVG.

```
# plot spatial expression of top-ranked SVG
ix <- which(rowData(spe)$rank == 1)
ix_name <- rowData(spe)$gene_name[ix]
ix_name

df <- as.data.frame(
  cbind(spatialCoords(spe), 
        expr = counts(spe)[ix, ]))

ggplot(df, aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres, 
               color = expr)) + 
  geom_point(size = 0.8) + 
  coord_fixed() + 
  scale_y_reverse() + 
  scale_color_gradient(low = "gray90", high = "blue", name = "counts") + 
  ggtitle(ix_name) + 
  theme_bw() + 
  theme(panel.grid = element_blank(), 
        axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank())
```


## Citation

Our paper describing `nnSVG` is available from bioRxiv:

- [Weber L.M. et al. (2022), "nnSVG: scalable identification of spatially variable genes using nearest-neighbor Gaussian processes", bioRxiv](https://www.biorxiv.org/content/10.1101/2022.05.16.492124v1)

