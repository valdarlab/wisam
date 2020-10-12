#' Weighted Genome Scan
#'
#' Performs a genome scan on heteroscedastic data.
#'
#' @param G A n by p matrix of genotypes, where n is the number of strains and p is the number of snps to be tested (can have missing values)
#' @param y A n length vector of mean phenotype for each strain
#' @param noise A n length vector of sample variances for each strain
#' @param counts A n length vector of number of replicates for each strain
#' @param X A n by q matrix of covariates (optional)
#' @param K A n by n genomic relationship matrix. Will be calculated if unspecified.
#' @param weights A string specifying the weights to be used. The following are permitted: "none", "samplevars", "limma", and "counts"
#'
#' @return A list containing:
#' \itemize{
#' \item{pvalue: A p length vector of p-values for every snp}
#' \item{ML0: A p length vector of log maximum likelihood under the null hypothesis}
#' \item{ML1: A p length vector of log maximum likelihood under the alternative hypothesis}
#' \item{beta: A p length vector of regression parameter estimate for slope under the alternative hypothesis}
#' \item{s2: A p length vector of the sum of variance components under the null hypothesis}
#' \item{h2: A p length vector of the heritability estimate under the null hypothesis}
#' }
#'
#' @import tidyverse
#' @import stringr
#' @importFrom limma squeezeVar
#' @importFrom emma emma.kinship
#'
#' @export
wisam <- function(G, y, noise = NULL, counts = NULL, X, K, weights = "none"){

  # number of strains
  n <- length(y)

  #### UNACCEPTABLE MISSINGNESS ####
  if (missing(y)) { stop('Must provide y (n-vector of phenotypes) to run a genome Scan.') }
  if (missing(G)) { stop('Must provide at least one snp to run a genome scan.')}

  #### ACCEPTABLE MISSINGNESS ####
  # initialize X to an intercept if missing
  if (missing(X)) { X <- matrix(data = 1, nrow = n) }
  # initialize K using the G matrix and emma package if missing
  if (missing(K)) { K <- emma.kinship(t(G), "additive", "all") }
  # intialize noise and counts if missing
  if (missing(noise)) { noise <- matrix(data = 1, nrow = n)}
  if (missing(counts)) { counts <- matrix(data = 1, nrow = n)}

  #### CONDITIONS THAT CAUSE AN ERROR ####
  if (!all(c(nrow(X), nrow(G), nrow(K)))) {
    stop("Input dimensions don't match.")
  }
  # checks length of phenotypes, sample variances, and counts, need one for each strain
  # if noise, counts left NULL, will stop
  if(!all(c(length(y), length(noise), length(counts)))){
    stop("Input dimensions don't match.")
  }

  ## check for strains with 0 variance and take them out
  if(weights %in% c("samplevars")){
    ind <- which(noise == 0| counts == 1)
    if (length(ind) > 0){
      y <- y[-ind]
      noise <- noise[-ind]
      counts <- counts[-ind]
      K <- K[-ind,-ind]
      G <- G[-ind,]
      X <- X[-ind,]
      print(ind)
    }
  }

  ######### WEIGHTS
  sample_vars = noise
  if (weights == "samplevars"){
    print("samplevars")
    weights = counts/sample_vars
    weights = weights/sum(weights)*sum(counts)
  } else if (weights == "limma"){
    print("limma")
    vars_shrink_limmar = squeezeVar(sample_vars, counts-1, robust = TRUE)$var.post
    weights = counts/vars_shrink_limmar
    weights = weights/sum(weights)*sum(counts)
  } else if(weights == "none"){
    print("no weights")
    weights <- rep(1, dim(K)[1])
  } else if(weights == "counts"){
    print("counts")
    weights <- counts
    weights = weights/sum(weights)*sum(counts)
  }

  ######## Find unique SNPs
  temp = apply(G, 2, mapping)
  strings = apply(temp, 2, paste0, collapse = "")
  uniques_indices = tapply(seq_along(strings), strings, identity)[unique(strings)]
  uniques = unique(strings)
  unique_counts = lapply(uniques_indices, length) %>% unname() %>% unlist()

  # run the scan (see functions.R)
  p_value_het = scan.strain.means(uniques, y, X, K, weights)
    ps = data.frame(indices = uniques_indices %>% unlist %>% unname(),
                    pvalue = rep(p_value_het$ps, unique_counts),
                    ML1 = rep(p_value_het$theta1, unique_counts),
                    ML0 = rep(p_value_het$theta0, unique_counts),
                    beta = rep(p_value_het$beta, unique_counts),
                    h2 = rep(p_value_het$h20, unique_counts),
                    s2 = rep(p_value_het$s20, unique_counts)) %>% dplyr::arrange(indices) %>%
      dplyr::select(-indices) #%>% unlist() %>% unname()
  ps %>% as.list()
}
