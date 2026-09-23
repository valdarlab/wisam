#' Posterior Mean Estimator of the Variance with Integrated Combined Likelihood
#'
#' Performs a genome scan on heteroscedastic data.
#'
#' @param G A s by p matrix of genotypes, where s is the number of strains and p is the number of snps to be tested. This matrix can have missing values, and each SNP should be coded as 0, 1, 0.5 or NA
#' @param y A N length vector of phenotypes for each individual organism, where N is the total number of individuals
#' @param strains A s by N incidence matrix that maps every individual to a strain
#' @param X A s by q matrix of covariates (optional)
#' @param K A s by s genomic relationship matrix. Will be calculated if unspecified.
#' @param weights A string specifying the weights to be used. The following are permitted: "none", "samplevars","EB", "limma", "counts", and "user"
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
combined_likelihood <- function(pheno_long,N,ns){
  combined_ss_list <- rep(NA,N)
  pheno_means = pheno_long %>% dplyr::group_by(strains) %>% dplyr::summarise(mean = mean(y),
                                                                             noise = var(y),
                                                                             counts = dplyr::n())
  for(i in 1:N){
    strain <- unique(pheno_long$strains)[i]
    mean_y <- pheno_means$mean[i]
    temp <- (pheno_long[pheno_long$strains == strain,1]- mean_y)^2
    combined_ss_list[i] <- sum(temp,na.rm = TRUE)
  }
  return(sum(combined_ss_list)/ sum(ns-1))
}
#'
#' @export
# front facing function to shrink variances
zshrink <- function(pheno_long, plot_densities = FALSE, plot_file = NA){
  pheno_summary <- pheno_long %>% dplyr::group_by(strains) %>% dplyr::summarise(mean = mean(y),
                                                                             noise = var(y),
                                                                             counts = dplyr::n())

  y <- pheno_summary$mean
  sigma2_observed <- pheno_summary$noise
  ns <- pheno_summary$counts
  N <- length(ns)
  # estimating hyperparameters
  sigma2_bar <- combined_likelihood(pheno_long,N,ns)
  #excluding NA (cases with one observation) from the estimates for lambda and sigmabar
  lambda <- var(sigma2_observed,na.rm = TRUE)
  a <- (sigma2_bar^2)/lambda + 2
  b <- sigma2_bar* (sigma2_bar^2/lambda + 1)

  # Posterior Mean Estimate
  sigma2_estimate <- rep(NA,N)
  for(i in 1:N){
    #catch cases with only 1 observation. otherwise ns = 1 and sigma2_observe = NA
    if(ns[i] == 1){
      numerator <- b
      denominator <- a - 1
    }else{
      numerator <- b + (ns[i]-1)* sigma2_observed[i] / 2
      denominator <- a + (ns[i]-1)/2 - 1
    }
    sigma2_estimate[i] <- numerator / denominator
  }

  # if(plot_densities){
  #   jpeg(filename = plot_file)
  #   den_observed <- density(sigma2_observed)
  #   den_estimate <- density(sigma2_estimate)
  #   xmin <- (min(min(den_estimate$x), min(den_observed)))
  #
  #   plot(density(den_observed),col = "black",xlim = c(),main = "Variance Shrinkage of Simulated Normal Data")
  #   lines(density(sigma2_observed))
  #   legend(x = 0.8, , legend = c("Observed", "Estimated"),
  #          fill = c("black", "red"),cex = .8)
  #
  # }
  return(data.frame("VarObserved" = sigma2_observed,"VarEstimated" = sigma2_estimate, "NIndividuals" = ns))
}
