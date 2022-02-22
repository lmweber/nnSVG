# nnSVG

[![R build status](https://github.com/lmweber/nnSVG/workflows/R-CMD-check-bioc/badge.svg)](https://github.com/lmweber/nnSVG/actions)

R/Bioconductor package implementing our method `nnSVG` for scalable identification of spatially variable genes (SVGs) using nearest-neighbor Gaussian processes in spatially resolved transcriptomics (ST) data.

The package is integrated into Bioconductor and will be submitted to Bioconductor soon.


## Installation

The development version of the package can be installed from GitHub:

```
remotes::install_github("lmweber/nnSVG")
```

The current development version also depends on the development version of the [BRISC](https://cran.r-project.org/package=BRISC) R package, which can be installed from the following repository:

```
remotes::install_github("ArkajyotiSaha/BRISC-extensions")
```


## Tutorial

An extended example and tutorial is available in the package vignette.


## Citation

A preprint describing `nnSVG` will be submitted to bioRxiv soon.


## Additional references

- Nearest-neighbor Gaussian processes (NNGP): [Datta et al. (2016)](https://www.tandfonline.com/doi/full/10.1080/01621459.2015.1044091)
- BRISC: [Saha and Datta (2018)](https://onlinelibrary.wiley.com/doi/full/10.1002/sta4.184)

