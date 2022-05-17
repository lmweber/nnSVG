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


## Citation

Our paper describing `nnSVG` is available from bioRxiv:

- [Weber L.M. et al. (2022), "nnSVG: scalable identification of spatially variable genes using nearest-neighbor Gaussian processes", bioRxiv](https://www.biorxiv.org/content/10.1101/2022.05.16.492124v1)



## Workflow

An example workflow is available in the package vignette from [Bioconductor](https://bioconductor.org/packages/nnSVG). A short overview is also provided below.


**Load packages**

```
library(nnSVG)
library(STexampleData)
library(scran)
library(ggplot2)
```


**Load example dataset**

```
# load example dataset from STexampleData package
spe <- Visium_humanDLPFC()
dim(spe)
```

```
## [1] 33538  4992
```


**Preprocessing**

```
# keep spots over tissue
spe <- spe[, colData(spe)$in_tissue == 1]
dim(spe)
```

```
## [1] 33538  3639
```

```
# spot-level quality control: already performed on this example dataset

# filter low-expressed and mitochondrial genes
# using function from nnSVG package with default filtering parameters
spe <- filter_genes(spe)
```

```
## Gene filtering: removing mitochondrial genes
## removed 13 mitochondrial genes
## Gene filtering: retaining genes with at least 3 counts in at least 0.5% (n = 19) of spatial locations
## removed 30216 out of 33525 genes due to low expression
```

```
# calculate log-transformed normalized counts (logcounts) using scran package
set.seed(123)
qclus <- quickCluster(spe)
spe <- computeSumFactors(spe, cluster = qclus)
spe <- logNormCounts(spe)
assayNames(spe)
```

```
## [1] "counts"    "logcounts"
```


**Subset data for this example**

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

```
## [1]   16 3639
```


**Run nnSVG**

```
# set seed for reproducibility
set.seed(123)
# run nnSVG using a single thread for this example workflow
spe <- nnSVG(spe, n_threads = 1)

# show results
rowData(spe)
```

```
## DataFrame with 16 rows and 17 columns
```


**Investigate results**

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
```

```
## 
## FALSE  TRUE 
##     7     9
```


```
# show results for top n SVGs
n <- 10
rowData(spe)[order(rowData(spe)$rank)[1:n], ]
```

```
## DataFrame with 10 rows and 17 columns
##                         gene_id   gene_name    feature_type   sigma.sq
##                     <character> <character>     <character>  <numeric>
## ENSG00000168314 ENSG00000168314        MOBP Gene Expression 1.86459294
## ENSG00000132639 ENSG00000132639      SNAP25 Gene Expression 0.33829228
## ENSG00000211592 ENSG00000211592        IGKC Gene Expression 0.59161928
## ENSG00000244734 ENSG00000244734         HBB Gene Expression 0.34819123
## ENSG00000183036 ENSG00000183036        PCP4 Gene Expression 0.22354847
## ENSG00000122585 ENSG00000122585         NPY Gene Expression 0.29511061
## ENSG00000129562 ENSG00000129562        DAD1 Gene Expression 0.03687246
## ENSG00000114923 ENSG00000114923      SLC4A3 Gene Expression 0.01123674
## ENSG00000133606 ENSG00000133606       MKRN1 Gene Expression 0.00543859
## ENSG00000149923 ENSG00000149923       PPP4C Gene Expression 0.12004347
##                    tau.sq        phi    loglik   runtime      mean       var
##                 <numeric>  <numeric> <numeric> <numeric> <numeric> <numeric>
## ENSG00000168314  0.371646   0.922937  -3716.46     0.959  0.841100  1.382681
## ENSG00000132639  0.440346   3.570016  -3940.04     0.668  3.464790  0.779762
## ENSG00000211592  0.464762  20.035566  -4580.98     1.618  0.630200  1.042847
## ENSG00000244734  0.365750  27.611193  -4114.59     2.233  0.418996  0.729640
## ENSG00000183036  0.456889   8.700988  -4041.98     1.001  0.684281  0.681316
## ENSG00000122585  0.302841  68.183198  -4087.69     1.375  0.401353  0.599801
## ENSG00000129562  0.484816   8.805056  -3942.80     2.051  0.561114  0.523034
## ENSG00000114923  0.237750  16.239042  -2621.78     0.992  0.249525  0.249055
## ENSG00000133606  0.277671   0.537947  -2861.92     1.873  0.295867  0.283165
## ENSG00000149923  0.132992 198.872410  -2660.25     5.132  0.235632  0.253096
##                     spcov   prop_sv loglik_lm    LR_stat      rank        pval
##                 <numeric> <numeric> <numeric>  <numeric> <numeric>   <numeric>
## ENSG00000168314  1.623470 0.8338075  -5752.58 4072.24582         1 0.00000e+00
## ENSG00000132639  0.167868 0.4344665  -4710.39 1540.68468         2 0.00000e+00
## ENSG00000211592  1.220514 0.5600433  -5239.35 1316.74199         3 0.00000e+00
## ENSG00000244734  1.408312 0.4877029  -4589.50  949.82287         4 0.00000e+00
## ENSG00000183036  0.690958 0.3285363  -4464.82  845.69334         5 0.00000e+00
## ENSG00000122585  1.353525 0.4935363  -4232.97  290.54324         6 0.00000e+00
## ENSG00000129562  0.342215 0.0706791  -3983.78   81.97837         7 0.00000e+00
## ENSG00000114923  0.424820 0.0451298  -2633.77   23.96612         8 6.24917e-06
## ENSG00000133606  0.249257 0.0192102  -2867.31   10.77036         9 4.58402e-03
## ENSG00000149923  1.470399 0.4744140  -2663.05    5.60524        10 6.06510e-02
##                        padj
##                   <numeric>
## ENSG00000168314 0.00000e+00
## ENSG00000132639 0.00000e+00
## ENSG00000211592 0.00000e+00
## ENSG00000244734 0.00000e+00
## ENSG00000183036 0.00000e+00
## ENSG00000122585 0.00000e+00
## ENSG00000129562 0.00000e+00
## ENSG00000114923 1.24983e-05
## ENSG00000133606 8.14937e-03
## ENSG00000149923 9.70416e-02
```


**Plot expression**

Plot expression of the top-ranked SVG.

```
# plot spatial expression of top-ranked SVG
ix <- which(rowData(spe)$rank == 1)
ix_name <- rowData(spe)$gene_name[ix]
ix_name
```

```
## [1] "MOBP"
```


```
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

<img src="expression_top_SVG.png" alt="Spatial expression plot of top-ranked SVG" title="Spatial expression plot of top-ranked SVG" width="400px">

