# wisam
Weighted Inbred Strain Association Mapping -- an R software package

When conducting Genome Wide Association Studies with inbred organisms, each strain may contribute an unequal amount of evironmental variance, violating the homogeneity of variances assumption in linear regression. The `wisam` package handles variance heterogeneity by weighting each strain by a standardized inverse standard error (where standard error refers to the sample variance / number of replicates). `wisam` is also capable of more stably estimating the sample variance in the case of small replicate sizes, by using an Empirical Bayes shrinkage method from the `limma` package in R. This can be implemented in `wisam` by choosing `weights = limma`.

## Installation 

```
install.packages("devtools")
devtools::install_github('williamvaldar/wisam')
```

## Usage
```
y = rnorm(1000)
G = matrix(rbinom(1000,1,0.5), nrow = 100)
strains = rep(1:100, each = 10) 
wisam::wisam(G, y, strains, weights = "limma")
```
