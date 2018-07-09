# R package: slingshot
[![Build Status](https://travis-ci.org/kstreet13/slingshot.svg?branch=master)](https://travis-ci.org/kstreet13/slingshot)
[![Coverage Status](https://img.shields.io/codecov/c/github/kstreet13/slingshot/master.svg)](https://codecov.io/github/kstreet13/slingshot?branch=master)

Provides functions for inferring continuous, branching lineage structures in low-dimensional data. Slingshot was designed to model developmental trajectories in single-cell RNA sequencing data and serve as a component in an analysis pipeline after dimensionality reduction and clustering. It is flexible enough to handle arbitrarily many branching events and allows for the incorporation of prior knowledge through supervised graph construction.

## Installation

```r
source("https://bioconductor.org/biocLite.R")
biocLite("kstreet13/slingshot")
```

## Issues and bug reports

Please use https://github.com/kstreet13/slingshot/issues to submit issues, bug reports, and comments.
