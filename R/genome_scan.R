#' Weighted Genome Scan
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
#' \item{ML0: A p length vector of log maximum likelihod under the null hypothesis}
#' \item{ML1: A p length vector of log maximum likelihood under the alternative hypothesis}
#' \item{beta: A p length vector of regression parameter estimate for slope under the alternative hypothesis}
#' \item{weights: A p length vector of regression weights estimated in wisam}
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
wisam <- function(G, y, strains, X, K, weights = "none",
                  user_weights = NULL, export_weights = FALSE, export_cooks= FALSE, export_s2hat = FALSE,use_individuals = FALSE,verbose = FALSE){
  # number of strains

  if(is.null(dim(strains))){
    s <- length(unique(strains))
    if(length(y) != length(strains) ){
      stop("Input dimensions don't match.")
    }
  }else{
    s <- nrow(strains)
    if(length(y) != ncol(strains) ){
      stop("Input dimensions don't match.")
    }
  }

  #### UNACCEPTABLE MISSINGNESS ####
  if (missing(y)) { stop('Must provide y (vector of phenotypes) to run a genome Scan.') }
  if (missing(G)) { stop('Must provide at least one snp to run a genome scan.')}
  if (is.null(user_weights) & weights == "user") { stop('Must provide user_weights (vector of weights) if weights = "user".')}
  if(!(weights %in% c("none", "samplevars", "zshrink", "limma","limma_nr","vashr","vashr_single", "counts", "user"))){stop('Weights must be one of "none", "samplevars","zshrink","limma_nr", "limma", "vashr", "vashr_single", "counts", and "user"')}
  if(export_s2hat & weights %in% c("none", "samplevars", "counts", "user")){stop("No variance shrinkage applied, cannot return s2hat.")}
  #### ACCEPTABLE MISSINGNESS ####
  # initialize K using the G matrix and emma package if missing
  if (missing(K)) { K <- emma.kinship(t(G), "additive", "all") }

  # initialize X to an intercept if missing
  if (missing(X) ) { X <- matrix(data = 1, nrow = nrow(K)) }
  #### CONDITIONS THAT CAUSE AN ERROR ####
  # if (!all(sapply(list(nrow(X), nrow(G), nrow(K)),
  #                 FUN = identical, nrow(strains)))){
  #   stop("Input dimensions don't match.")
  # }
  # checks length of phenotypes and strains

  # if(nrow(G) != nrow(strains)){
  #   stop("Input dimensions don't match.")
  # }
  # check that G is coded correctly
  G_unique <- as.vector(as.matrix(G)) %>% unique()
  if(!setequal(G_unique, c(0,1,NA,0.5)) & !setequal(G_unique, c(0,1,NA)) &
     !setequal(G_unique, c(0,1,0.5)) & !setequal(G_unique, c(0,1))){
    #stop("Each SNP should be coded as 0, 1, 0.5, or NA")
  }
  #check for order of y
  if(is.null(names(y))){
    warning("y vector does not have names - please be sure that the y vector is in the same order as K and G.")
  }
  if(!is.null(dim(strains))){ strains = apply(strains, 2, function(x) x*c(1:nrow(strains))) %>% colSums()
  }


  pheno_long = data.frame(y = y, strains = strains)
  pheno_means = pheno_long %>% dplyr::group_by(strains) %>% dplyr::summarise(mean = mean(y),
                                                               noise = var(y),
                                                               counts = dplyr::n())
  pheno_means <- as.data.frame(pheno_means)
  pheno_long$noise = pheno_means[match(pheno_long$strains,pheno_means$strains),"noise"]
  pheno_long$counts = pheno_means[match(pheno_long$strains,pheno_means$strains),"counts"]
  if(use_individuals){
    y = pheno_long$y
    noise = pheno_long$noise
    counts = pheno_long$counts
  }else{
    y = pheno_means$mean
    noise = pheno_means$noise
    counts = pheno_means$counts
    names = pheno_means$strains
  }


  ## check for strains with 0 variance and take them out for sample variance estimates
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
  switch(EXPR = weights,
         samplevars = {
           if(verbose) print("caluclating samplevars")
           weights = mean(sample_vars)/sample_vars
           ds <-  1/(weights * counts)
           },
         zshrink = {
           if(verbose) print("calculating shrinkage estimates with integrated conditional likelihood ebayes")
           shrink_estimates <- zshrink(pheno_long, strains)$VarEstimated
           if(use_individuals) shrink_estimates = rep(shrink_estimates, pheno_means$counts)
           weights = (mean(shrink_estimates) / shrink_estimates )
           ds <-  1/(weights * counts)
           #weights = weights/(sum(weights)*sum(counts))
         },limma_nr = {
           if(verbose) print("calculating shrinkage estimates with non-robust limma")
           shrink_estimates = squeezeVar(sample_vars, counts-1, robust = FALSE)$var.post
           weights = mean(shrink_estimates) / shrink_estimates
           ds <- 1/(weights * counts)
         },
         limma = {
           if(verbose) print("calculating shrinkage estimates with robust limma")
           shrink_estimates = squeezeVar(sample_vars, counts-1, robust = TRUE)$var.post
           weights = (mean(shrink_estimates) / shrink_estimates )
           ds <- 1/(weights * counts)
         },
         vashr = {
           if(verbose) print("calculating shrinkage estimates with vashr")
           shrink_estimates = vashr::vash(sample_vars, df = counts[1]-1)$sd.post
           weights = mean(shrink_estimates) / shrink_estimates
           ds <- 1/(weights * counts)
         },
         vashr_single = {
           #vashr cannot handle NA
           sample_vars[is.na(sample_vars)] <- 0
           if(verbose) print("calculating shrinkage estimates with vashr and one prior")
           shrink_estimates = vashr::vash(sample_vars,singlecomp = TRUE, df = counts[1]-1)$sd.post
           weights = mean(shrink_estimates) / shrink_estimates
           ds <- 1/(weights * counts)
         },
         none = {
           if(verbose) print("no weights used")
           weights <- rep(1, dim(K)[1])
           ds <- weights
         },
         counts = {
           if(verbose) print("using counts as weights")
           weights <- counts
           ds <- 1/(weights)
         },
         user = {
           weights <- user_weights
           ds <- 1/(weights * counts)
         },
         stop("unkown weights input")
         )

#
#   if (weights == "samplevars"){
#     print("samplevars")
#     weights = counts/sample_vars
#     weights = weights/sum(weights)*sum(counts)
#   } else if (weights == "eb"){
#   } else if (weights == "limma"){
#
#   } else if(weights == "none"){
#
#   } else if(weights == "counts"){
#
#   } else if(weights == "user"){
#
#   }

  ######## Find unique SNPs
  temp = apply(G, 2, mapping)
  strings = apply(temp, 2, paste0, collapse = "")
  uniques_indices = tapply(seq_along(strings), strings, identity)[unique(strings)]
  uniques = unique(strings)
  G_names_uniques <- colnames(G)[unname(unlist(uniques_indices))]
  unique_counts = lapply(uniques_indices, length) %>% unname() %>% unlist()

  if(export_cooks){
    results = scan.strain.means.weight.cookd(G = uniques, y = y, X = X, K = K, weights = ds)
    print("successfully returned from scan strain means")
    cooksD = results$cooksD
    resids = results$results
    weights = results$weights
    results = results$results
  }else{
    results = scan.strain.means(uniques, y, X, K, ds)
  }

  if(!export_s2hat){
    s2return <- sample_vars
  }else{
    s2return <- shrink_estimates
  }

   ps = data.frame(indices = uniques_indices %>% unlist %>% unname(),
                  pvalue = rep(results$ps, unique_counts),
                  ML1 = rep(results$theta1, unique_counts),
                  ML0 = rep(results$theta0, unique_counts),
                  beta = rep(results$beta, unique_counts),
                  h2 = rep(results$h20, unique_counts),
                  s2 = rep(results$s20, unique_counts)) %>% dplyr::arrange(indices) %>%
    dplyr::select(-indices) #%>% unlist() %>% unname()
    rownames(ps) <- G_names_uniques

    if(export_weights & export_cooks){
      return(list("results" = ps, "weights" = weights, "s2s" = s2return,
                  "cooksdistance" = cooksD, "residuals" = resids))
    }else if(export_cooks){
      return(list("results" = ps,"cooksdistance" = cooksD, "s2s" = s2return))
    }else if(export_weights ){
      return(list("results" = ps, "weights" = weights, "s2s" = s2return))
    }

    list(ps, "s2s" = s2return)

}
