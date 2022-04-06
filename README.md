# nnSVG

[![R build status](https://github.com/lmweber/nnSVG/workflows/R-CMD-check-bioc/badge.svg)](https://github.com/lmweber/nnSVG/actions)

R/Bioconductor package implementing our method `nnSVG` for scalable identification of spatially variable genes (SVGs) in spatially resolved transcriptomics (ST) data.

`nnSVG` is based on nearest-neighbor Gaussian processes ([Datta et al., 2016](https://www.tandfonline.com/doi/full/10.1080/01621459.2015.1044091), [Finley et al., 2019](https://www.tandfonline.com/doi/full/10.1080/10618600.2018.1537924)) and uses the BRISC algorithm ([Saha and Datta, 2018](https://onlinelibrary.wiley.com/doi/full/10.1002/sta4.184)) for model fitting and parameter estimation. `nnSVG` allows identification and ranking of SVGs with flexible length scales across a tissue slide or within spatial domains defined by covariates. The method scales linearly with the number of spatial locations and can be applied to datasets containing thousands or more spatial locations.

The `nnSVG` package is integrated into the Bioconductor framework and uses the [SpatialExperiment](https://bioconductor.org/packages/SpatialExperiment) class to store ST data and results.

More details on the method will be provided in our paper (to be submitted to bioRxiv soon).

The package is currently available from GitHub and has been submitted to Bioconductor.


## Installation

The development version of the package can be installed from GitHub as follows, including updated versions of dependencies. (Note the `ref = "release"` argument.)

```
remotes::install_github("drighelli/SpatialExperiment")
remotes::install_github("lmweber/STexampleData")
install.packages("BRISC")
remotes::install_github("lmweber/nnSVG", ref = "release")
```


## Tutorial

An extended example and tutorial is available in the package vignette.


## Citation

A preprint describing `nnSVG` will be submitted to bioRxiv soon.


## Additional references

- Nearest-neighbor Gaussian processes (NNGP): [Datta et al. (2016)](https://www.tandfonline.com/doi/full/10.1080/01621459.2015.1044091)
- BRISC: [Saha and Datta (2018)](https://onlinelibrary.wiley.com/doi/full/10.1002/sta4.184)

