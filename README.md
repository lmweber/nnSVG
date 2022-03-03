# nnSVG

[![R build status](https://github.com/lmweber/nnSVG/workflows/R-CMD-check-bioc/badge.svg)](https://github.com/lmweber/nnSVG/actions)

R/Bioconductor package implementing our method `nnSVG` for scalable identification of spatially variable genes (SVGs) using nearest-neighbor Gaussian processes in spatially resolved transcriptomics (ST) data.

The package is integrated into the Bioconductor framework and uses the [SpatialExperiment](https://bioconductor.org/packages/SpatialExperiment) class to store ST data and results.

The package will be submitted to Bioconductor soon.


## Installation

The development version of the package can be installed from GitHub:

```
remotes::install_github("lmweber/nnSVG")
```


## Tutorial

An extended example and tutorial is available in the package vignette.


## Citation

A preprint describing `nnSVG` will be submitted to bioRxiv soon.


## Additional references

- Nearest-neighbor Gaussian processes (NNGP): [Datta et al. (2016)](https://www.tandfonline.com/doi/full/10.1080/01621459.2015.1044091)
- BRISC: [Saha and Datta (2018)](https://onlinelibrary.wiley.com/doi/full/10.1002/sta4.184)

