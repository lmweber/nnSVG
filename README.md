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


## Tutorial

An extended example and tutorial is available in the package vignette.


## Citation

Our paper describing `nnSVG` is available from bioRxiv:

- [Weber L.M. et al. (2022), "nnSVG: scalable identification of spatially variable genes using nearest-neighbor Gaussian processes", bioRxiv](https://www.biorxiv.org/content/10.1101/2022.05.16.492124v1)

