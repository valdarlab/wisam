# wisam
Weighted Inbred Strain Association Mapping -- an R software package

When conducting Genome Wide Association Studies with inbred organisms, each strain may contribute an unequal amount of environmental variance, violating the homogeneity of variances assumption in linear regression. The `wisam` package handles variance heterogeneity by weighting each strain by a standardized inverse standard error (where standard error refers to the sample variance / number of replicates). `wisam` is also capable of more stably estimating the sample variance in the case of small replicate sizes, by using an Empirical Bayes shrinkage method from the `limma` package in R. This can be implemented in `wisam` by choosing `weights = "limma"`.

## Installation 

Use the following code chunk to install `wisam`.
```
if(!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("limma")

install.packages("devtools")
devtools::install_github('valdarlab/wisam')
```

## Usage
The following code chunk generates a toy example dataset to show how the function is used.

```
library(wisam)
library(tidyverse)
y = rnorm(1000)
G = matrix(rbinom(1000,1,0.5), nrow = 100)
strains = rep(1:100, each = 10)
tab <- table(cbind.data.frame(ID = 1:length(strains), 
                              Strain = strains)) %>% as.data.frame()
Z = tidyr::spread(tab, Strain, Freq)[-1] %>% t()
wisam::wisam(G, y, Z, weights = "limma")
```
