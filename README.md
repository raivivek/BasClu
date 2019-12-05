## BasClu

R package for "Hierarchical Bayesian Model for Clustering of Single-cell RNA-seq Data" as published in Liu et al (2019).

## Installation
This package can be installed from github. It requires "Rcpp" and "RcppArmadillo" as
dependencies, both of which can are installed as dependencies.

```
devtools::install_github("raivivek/BasClu")
```

This is a work in progress. Use at your own risk.

## Example
```
library(BasClu)

data("okaty")
dim(okaty)

idx_genes <- sample(1:nrow(okaty), 100)

count_mat <- as.matrix(okaty[, -(1:3)])
transcript_lengths <- matrix(okaty[, 2], nrow = 1)
result <- BasClu(count_mat, transcript_lengths, clustering = "BasCluZ", ...)
```

### Example dataset
Package includes scRNA-seq example data from Okaty et al. (2015) consisting of a
transcriptomic profiling of 5HT neurons.

## Reference
Liu, Y., Warren, J. L. & Zhao, H. A hierarchical Bayesian model for single-cell
clustering using RNA-sequencing data. Ann. Appl. Stat. 13, 1733–1752 (2019).
