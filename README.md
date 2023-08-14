# nnSVG

[![R build status](https://github.com/lmweber/nnSVG/workflows/R-CMD-check-bioc/badge.svg)](https://github.com/lmweber/nnSVG/actions)


## Overview

`nnSVG` is a method for scalable identification of spatially variable genes (SVGs) in spatially-resolved transcriptomics data.

The `nnSVG` method is based on nearest-neighbor Gaussian processes ([Datta et al., 2016](https://www.tandfonline.com/doi/full/10.1080/01621459.2015.1044091), [Finley et al., 2019](https://www.tandfonline.com/doi/full/10.1080/10618600.2018.1537924)) and uses the BRISC algorithm ([Saha and Datta, 2018](https://onlinelibrary.wiley.com/doi/full/10.1002/sta4.184)) for model fitting and parameter estimation. `nnSVG` allows identification and ranking of SVGs with flexible length scales across a tissue slide or within spatial domains defined by covariates. The method scales linearly with the number of spatial locations and can be applied to datasets containing thousands or more spatial locations.

`nnSVG` is implemented as an R package within the Bioconductor framework, and is available from [Bioconductor](https://bioconductor.org/packages/nnSVG).

Our paper describing the method is available from [Nature Communications](https://www.nature.com/articles/s41467-023-39748-z).


## Installation

The package can be installed from [Bioconductor](https://bioconductor.org/packages/nnSVG) as follows, using R version 4.2 or above:

```r
install.packages("BiocManager")
BiocManager::install("nnSVG")
```

Alternatively, the latest development version of the package can also be installed from GitHub:

```r
remotes::install_github("lmweber/nnSVG")
```

If you are installing from GitHub, the following dependency packages may need to be installed manually from Bioconductor and CRAN (these are installed automatically if you install from Bioconductor instead):

```r
install.packages("BiocManager")
BiocManager::install("SpatialExperiment")
BiocManager::install("STexampleData")
install.packages("BRISC")
```


## Tutorial

A detailed tutorial is available in the package vignette from [Bioconductor](https://bioconductor.org/packages/nnSVG). A direct link to the tutorial / package vignette is available [here](https://bioconductor.org/packages/release/bioc/vignettes/nnSVG/inst/doc/nnSVG.html).


## Input data format

In the examples below, we assume the input data are provided as a [SpatialExperiment](https://bioconductor.org/packages/SpatialExperiment) Bioconductor object. In this case, the outputs are stored in the `rowData` of the `SpatialExperiment` object.

Alternatively, the inputs can also be provided as a numeric matrix of normalized and transformed counts (e.g. log-transformed normalized counts, also known as logcounts) and a numeric matrix of spatial coordinates.


## Example workflow

A short example workflow is shown below. This is a modified version of the full tutorial available in the package vignette from [Bioconductor](https://bioconductor.org/packages/nnSVG). A direct link to the tutorial / package vignette is available [here](https://bioconductor.org/packages/release/bioc/vignettes/nnSVG/inst/doc/nnSVG.html)).


**Load packages**

```r
library(nnSVG)
library(STexampleData)
library(scran)
library(ggplot2)
```


**Load example dataset**

```r
# load example dataset from STexampleData package
spe <- Visium_humanDLPFC()
dim(spe)
```

```r
## [1] 33538  4992
```


**Preprocessing**

```r
# keep spots over tissue
spe <- spe[, colData(spe)$in_tissue == 1]
dim(spe)
```

```r
## [1] 33538  3639
```

```r
# spot-level quality control: already performed on this example dataset
```

```r
# filter low-expressed and mitochondrial genes
# using function from nnSVG package with default filtering parameters
spe <- filter_genes(spe)
```

```r
## Gene filtering: removing mitochondrial genes
## removed 13 mitochondrial genes
## Gene filtering: retaining genes with at least 3 counts in at least 0.5% (n = 19) of spatial locations
## removed 30216 out of 33525 genes due to low expression
```

```r
# calculate logcounts (log-transformed normalized counts) using scran package
# using library size factors
spe <- computeLibraryFactors(spe)
spe <- logNormCounts(spe)
assayNames(spe)
```

```r
## [1] "counts"    "logcounts"
```


**Subset data for this example**

```r
# select small set of random genes and several known SVGs for faster runtime in this example workflow
set.seed(123)
ix_random <- sample(seq_len(nrow(spe)), 10)
known_genes <- c("MOBP", "PCP4", "SNAP25", "HBB", "IGKC", "NPY")
ix_known <- which(rowData(spe)$gene_name %in% known_genes)
ix <- c(ix_known, ix_random)

spe <- spe[ix, ]
dim(spe)
```

```r
## [1]   16 3639
```


**Run nnSVG**

```r
# set seed for reproducibility
# run nnSVG using a single thread for this example workflow
set.seed(123)
spe <- nnSVG(spe, n_threads = 1)

# show results
rowData(spe)
```

```r
## DataFrame with 16 rows and 17 columns
## [...]
```


**Investigate results**

The results are stored in the `rowData` of the `SpatialExperiment` object.

The main results of interest are:

- `LR_stat`: likelihood ratio (LR) statistics used to rank SVGs
- `rank`: rank of top SVGs according to LR statistics
- `pval`: approximate p-values
- `padj`: approximate p-values adjusted for multiple testing
- `prop_sv`: effect size defined as proportion of spatial variance

```r
# number of significant SVGs
table(rowData(spe)$padj <= 0.05)
```

```r
## 
## FALSE  TRUE 
##     7     9
```


```r
# show results for top n SVGs
n <- 10
rowData(spe)[order(rowData(spe)$rank)[1:n], ]
```

```r
## DataFrame with 10 rows and 17 columns
##                         gene_id   gene_name    feature_type   sigma.sq    tau.sq
##                     <character> <character>     <character>  <numeric> <numeric>
## ENSG00000168314 ENSG00000168314        MOBP Gene Expression 1.38739383  0.364188
## ENSG00000132639 ENSG00000132639      SNAP25 Gene Expression 0.43003959  0.430106
## ENSG00000211592 ENSG00000211592        IGKC Gene Expression 0.56564845  0.455042
## ENSG00000244734 ENSG00000244734         HBB Gene Expression 0.32942113  0.353754
## ENSG00000183036 ENSG00000183036        PCP4 Gene Expression 0.23102220  0.452735
## ENSG00000122585 ENSG00000122585         NPY Gene Expression 0.28567359  0.280173
## ENSG00000129562 ENSG00000129562        DAD1 Gene Expression 0.02389607  0.464723
## ENSG00000114923 ENSG00000114923      SLC4A3 Gene Expression 0.01147170  0.237260
## ENSG00000133606 ENSG00000133606       MKRN1 Gene Expression 0.00632248  0.272432
## ENSG00000143543 ENSG00000143543         JTB Gene Expression 0.07541566  0.463623
##                        phi    loglik   runtime      mean       var     spcov
##                  <numeric> <numeric> <numeric> <numeric> <numeric> <numeric>
## ENSG00000168314   1.102018  -3663.60     0.631  0.805525  1.205673  1.462248
## ENSG00000132639   3.033847  -3912.70     0.450  3.451926  0.857922  0.189973
## ENSG00000211592  20.107022  -4531.64     1.054  0.622937  1.007454  1.207340
## ENSG00000244734  27.814098  -4044.96     1.559  0.411262  0.697673  1.395587
## ENSG00000183036   8.272278  -4026.22     0.419  0.687961  0.684598  0.698656
## ENSG00000122585  71.653290  -3995.23     0.843  0.393975  0.567383  1.356646
## ENSG00000129562  10.141894  -3842.24     0.590  0.549318  0.489167  0.281410
## ENSG00000114923  12.765645  -2617.36     0.658  0.250768  0.248816  0.427112
## ENSG00000133606   0.082764  -2831.51     0.612  0.295404  0.278806  0.269171
## ENSG00000143543 119.721419  -4036.28     0.731  0.654919  0.539172  0.419318
##                   prop_sv loglik_lm    LR_stat      rank        pval        padj
##                 <numeric> <numeric>  <numeric> <numeric>   <numeric>   <numeric>
## ENSG00000168314 0.7920804  -5503.33 3679.46397         1 0.00000e+00 0.00000e+00
## ENSG00000132639 0.4999614  -4884.19 1942.98556         2 0.00000e+00 0.00000e+00
## ENSG00000211592 0.5541822  -5176.53 1289.77508         3 0.00000e+00 0.00000e+00
## ENSG00000244734 0.4821910  -4507.99  926.04573         4 0.00000e+00 0.00000e+00
## ENSG00000183036 0.3378716  -4473.57  894.68884         5 0.00000e+00 0.00000e+00
## ENSG00000122585 0.5048609  -4131.87  273.27818         6 0.00000e+00 0.00000e+00
## ENSG00000129562 0.0489053  -3861.98   39.49098         7 2.65854e-09 6.07667e-09
## ENSG00000114923 0.0461207  -2632.02   29.31376         8 4.31119e-07 8.62238e-07
## ENSG00000133606 0.0226812  -2839.08   15.15227         9 5.12539e-04 9.11181e-04
## ENSG00000143543 0.1399077  -4039.07    5.59664        10 6.09124e-02 9.74599e-02
```


**Plot expression of top SVG**

Plot expression of the top-ranked SVG.

```r
# plot spatial expression of top-ranked SVG
ix <- which(rowData(spe)$rank == 1)
ix_name <- rowData(spe)$gene_name[ix]
ix_name
```

```r
## [1] "MOBP"
```


```r
df <- as.data.frame(cbind(spatialCoords(spe), expr = counts(spe)[ix, ]))

ggplot(df, aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres, 
               color = expr)) + 
  geom_point(size = 0.8) + 
  coord_fixed() + 
  scale_y_reverse() + 
  scale_color_gradient(low = "gray90", high = "blue", 
                       trans = "sqrt", breaks = range(df$expr), 
                       name = "counts") + 
  ggtitle(ix_name) + 
  theme_bw() + 
  theme(plot.title = element_text(face = "italic"), 
        panel.grid = element_blank(), 
        axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank())
```

<img src="https://user-images.githubusercontent.com/8062417/179638201-a6d0cc21-a625-4899-8ab4-9b082b1d3a8c.png" alt="Spatial expression plot of top-ranked SVG" title="Spatial expression plot of top-ranked SVG" width="350px">


## Citation

Our paper describing `nnSVG` is available from Nature Communications:

- [Weber L.M. et al. (2023), "nnSVG for the scalable identification of spatially variable genes using nearest-neighbor Gaussian processes", Nature Communications](https://www.nature.com/articles/s41467-023-39748-z)

