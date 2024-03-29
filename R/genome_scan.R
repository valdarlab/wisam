#' Weighted Genome Scan
#'
#' Performs a genome scan on heteroscedastic data.
#'
#' @param G A s by p matrix of genotypes, where s is the number of strains and p is the number of snps to be tested. This matrix can have missing values, and each SNP should be coded as 0, 1, 0.5 or NA
#' @param y A N length vector of phenotypes for each individual organism, where N is the total number of individuals
#' @param strains A s by N incidence matrix that maps every individual to a strain
#' @param X A s by q matrix of covariates (optional)
#' @param K A s by s genomic relationship matrix. Will be calculated if unspecified.
#' @param weights A string specifying the weights to be used. The following are permitted: "none", "samplevars", "limma", "counts", and "user"
#' @param user_weights A s length vector of weights for each strain, used if weights = "user"
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
#' @import statmod
#' @import plyr
#' @importFrom limma squeezeVar
#'
#' @export
wisam <- function(G, y, strains, X, K, weights = "none", user_weights = NULL){

  # number of strains
  s <- nrow(strains)

  #### UNACCEPTABLE MISSINGNESS ####
  if (missing(y)) { stop('Must provide y (vector of phenotypes) to run a genome Scan.') }
  if (missing(G)) { stop('Must provide at least one snp to run a genome scan.')}
  if (is.null(user_weights) & weights == "user") { stop('Must provide user_weights (vector of weights) if weights = "user".')}
  if(!(weights %in% c("none", "samplevars", "limma", "counts", "user"))){stop('Weights must be one of "none", "samplevars", "limma", "counts", and "user"')}

  #### ACCEPTABLE MISSINGNESS ####
  # initialize X to an intercept if missing
  if (missing(X)) { X <- matrix(data = 1, nrow = s) }
  # initialize K using the G matrix and emma package if missing
  if (missing(K)) { K <- emma.kinship(t(G), "additive", "all") }

  #### CONDITIONS THAT CAUSE AN ERROR ####
  if (!all(sapply(list(nrow(X), nrow(G), nrow(K)),
                  FUN = identical, nrow(strains)))){
    stop("Input dimensions don't match.")
  }
  # checks length of phenotypes and strains
  if(length(y) != ncol(strains)){
    stop("Input dimensions don't match.")
  }
  if(nrow(G) != nrow(strains)){
    stop("Input dimensions don't match.")
  }
  # check that G is coded correctly
  G_unique <- as.vector(as.matrix(G)) %>% unique()
  if(!setequal(G_unique, c(0,1,NA,0.5)) & !setequal(G_unique, c(0,1,NA)) &
     !setequal(G_unique, c(0,1,0.5)) & !setequal(G_unique, c(0,1))){
    stop("Each SNP should be coded as 0, 1, 0.5, or NA")
  }

  strains = apply(strains, 2, function(x) x*c(1:nrow(strains))) %>% colSums() %>% unname()
  pheno_long = data.frame(y = y, strains = strains)
  pheno_means = pheno_long %>% dplyr::group_by(strains) %>% dplyr::summarise(mean = mean(y),
                                                               noise = var(y),
                                                               counts = dplyr::n())
  y = pheno_means$mean
  noise = pheno_means$noise
  counts = pheno_means$counts

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
  } else if(weights == "user"){
    weights <- user_weights
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
