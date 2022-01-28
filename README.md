# nnSVG

[![R build status](https://github.com/lmweber/nnSVG/workflows/R-CMD-check-bioc/badge.svg)](https://github.com/lmweber/nnSVG/actions)

R/Bioconductor package implementing our method `nnSVG` for scalable identification of spatially variable genes (SVGs) using nearest-neighbor Gaussian processes in spatially resolved transcriptomics (ST) data.


## Installation

The package can be installed from GitHub as follows.

```
remotes::install_github("lmweber/nnSVG")
```

Note this currently requires the development version of the BRISC package containing additional updates.

```
remotes::install_github("ArkajyotiSaha/BRISC-extensions")
```


## Tutorial

An extended example and tutorial is available in the package vignette.


## References

- Nearest neighbor Gaussian processes (NNGP): [Datta et al. (2016)](https://www.tandfonline.com/doi/full/10.1080/01621459.2015.1044091)
- BRISC: [Saha and Datta (2018)](https://onlinelibrary.wiley.com/doi/full/10.1002/sta4.184)

