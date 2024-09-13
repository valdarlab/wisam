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
combined_likelihood <- function(fulldata,N,ns){
  combined_ss_list <- rep(NA,N)
  for(i in 1:N){
    mean_y <- mean(fulldata[i,],na.rm = TRUE)
    temp <- unlist(lapply(fulldata[i,], FUN = function(x){(x - mean_y)^2}))
    combined_ss_list[i] <- sum(temp,na.rm = TRUE)
  }
  return(sum(combined_ss_list)/ sum(ns-1))
}
#'
#' @export
# front facing function to shrink variances
estimateVar <- function(fulldata, strains){
  y <- unlist(apply(fulldata[,-ncol(fulldata)],MARGIN = 2,FUN = mean))
  sigma2_observed <- unlist(apply(fulldata[,-ncol(fulldata)],MARGIN = 2,FUN = var))
  ns <- unlist(apply(fulldata[,-ncol(fulldata)],MARGIN = 2,FUN = function(x){return(length(na.omit(x)))}))
  N <- length(pheno_means$mean)
  # estimating hyperparameters
  sigma2_bar <- combined_likelihood(fulldata[,-ncol(fulldata)],N,ns)
  #excluding NA (cases with one observation) from the estimates for lambda and sigmabar
  lambda <- var(sigma2_observed,na.rm = TRUE)
  a <- (sigma2_bar^2)/lambda + 2
  b <- sigma2_bar* (sigma2_bar^2/lambda + 1)

  # Posterior Mean Estimate
  sigma2_estimate <- rep(NA,N)
  for(i in 1:N){
    numerator <- b + (ns[i]-1)* sigma2_observed[i] / 2
    denominator <- a + (ns[i]-1)/2 - 1
    sigma2_estimate[i] <- numerator / denominator
  }
  # if there are NA, we can assume these are 0. When sigma2_hat is zero, the estimate is sigma_bar
  sigma2_estimate[is.na(sigma2_estimate)] <- sigma2_bar
  return(data.frame("VarObserved" = sigma2_observed,"VarEstimated" = sigma2_estimate, "NIndividuals" = ns))
}
