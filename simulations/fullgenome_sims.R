### This script runs full genome simulations 
### was run on Longleaf (computing cluster at UNC) for each different phenotype and method 

## things i need in the directory 
# bash script 
# genotype info 
# wisam directory with functions 

### Command Line Arguments 
args <- commandArgs(TRUE)

h2 <- as.numeric(args[1]) # heritability constant 
effect_size <- as.numeric(args[2]) # effect size (was 0.1 here)
genome_file <- as.character(args[3]) # file of genome used
output_file <- as.character(args[6])
n1 <- as.integer(args[4]) # replicate number if constant, or lower bound if nonconstant
n2 <- as.integer(args[5]) # 0 if replicate number constant, or upper bound if nonconstant 
# etc etc 

########################### Libraries needed ############################

# module load r/3.4.1
# R CMD INSTALL emma_1.1.2.tar
# R CMD INSTALL SimHaploidPop_0.1.tar 

library(tidyverse)
library(emma)
library(limma)
library(wisam)
# source("wisam/functions.R")
# source("wisam/genome_scan.R")

### Load R Data function 
loadRData <- function(filename) {
  load(filename)
  get(ls()[ls() != "filename"])
}

### load genotype
G = loadRData(genome_file)

# dupilicity of SNP configurations (want the true SNP effect to be in a unique SNP)
hmdp_t <- G %>% as.matrix() %>% t() %>% as.data.frame()
d1 <- duplicated.data.frame(hmdp_t, fromLast = FALSE)
d2 <- duplicated.data.frame(hmdp_t, fromLast = TRUE)
# uniques <- which(d1|d2)

missing <- apply(G, 2, function(x) sum(is.na(x)))

uniques <- which((d1|d2) & missing == 0)

### simulate phenotype 
  num_strain = nrow(G)
  variance = sample(2^c(-2, -1, 0, 1, 2), num_strain, replace = TRUE)
  variance_hom = rep(1, num_strain)
  if(n2 != 0) counts = sample(n1:n2, num_strain, replace = TRUE) else counts = rep(n1, num_strain)
  counts.df <- data.frame(strain  = 1:num_strain, count = counts)
  
  ### variance components 
  tau2 = h2 # to confirm that tau2 + sigma2 = 1 
  sig2 = tau2*(1-h2)/h2 

  ### Kinship Matrix via EMMA
  # dist <- dist(G, method = "manhattan") %>% as.matrix()
  # K <- 1 - dist/(2*ncol(G))
  K <- emma::emma.kinship(t(G), "additive", "all")
  
  ### Random Effect Matrix 
  A <- data.frame(value = MASS::mvrnorm(n = 1,
                                      mu = rep(0, num_strain),
                                      Sigma = K*tau2), count = counts.df$count)
A_large <- apply(A, 1, function(x) rep(x[1], x[2]))# %>% unlist() %>% unname()

### Random Error Matrices
E <- mapply(rnorm, n = counts, mean = 0, 
            sd = sqrt(variance*sig2)) 
E_large <- E #%>% unlist()

E_hom <- mapply(rnorm, n = counts, mean = 0, 
            sd = sqrt(variance_hom*sig2)) 
E_large_hom <- E_hom #%>% unlist()

if(length(unique(counts)) == 1){
  A_large <- A_large %>% as.data.frame() %>% as.list()
  E_large <- E_large %>% as.data.frame() %>% as.list()
  E_large_hom <- E_large_hom %>% as.data.frame() %>% as.list()
}


### Creation of incidence matrix 
strains <- apply(counts.df, 1, function(x) rep(x[1], x[2])) %>% as.list() %>% 
  unlist() %>% unname()   
tab <- table(cbind.data.frame(ID = 1:length(strains), 
                              Strain = strains)) %>% as.data.frame()
Z <- spread(tab, Strain, Freq)[-1] %>% as.matrix()


SNP = sample(uniques, 1) # pick a random SNP from the non-duplicated ones 
### SNP effect size function 
G_es = Z%*%(G[,SNP]%>%as.matrix())
A_es = A_large%>%unlist()
E_es = E_large %>% unlist()
E_es_hom = E_large_hom %>% unlist()

beta_het = sqrt(effect_size*sum((A_es + E_es - mean(A_es + E_es))^2)/
                  (1-effect_size)/sum((G_es-mean(G_es))^2))
beta_hom = sqrt(effect_size*sum((A_es + E_es - mean(A_es + E_es))^2)/
                  (1-effect_size)/sum((G_es-mean(G_es))^2))

### Phenotypes 
YL <- mapply(function(x, y, z) beta_het*z+ x + y,A_large,E_large, G[,SNP],
             SIMPLIFY = TRUE)
YL_hom <- mapply(function(x, y, z) beta_hom*z+ x + y,A_large,E_large_hom, G[,SNP],
             SIMPLIFY = TRUE)

if(length(unique(counts)) == 1){
  YL <- YL %>% as.data.frame() %>% as.list()
  YL_hom <- YL_hom %>% as.data.frame() %>% as.list()
}

### For strain means 
noise = lapply(YL, var) %>% unlist() %>% unname()
Y <- lapply(YL, mean) %>% unlist() %>% unname()
noise_hom = lapply(YL_hom, var) %>% unlist() %>% unname()
Y_hom <- lapply(YL_hom, mean) %>% unlist() %>% unname()

YL <- YL %>% unlist()
YL_hom <- YL_hom %>% unlist()

print("phenotypes generated")

### run scans 

############################## VARIANCE SHRINKAGE ################################
## this happens in the wisam package if selected 

## variance shrinkage 
#vars_shrink_limmar = squeezeVar(noise, counts, robust = TRUE)$var.post
#known_vars <- variance

#sample_vars_hom <- lapply(obs_phenotypes_tot, var) %>% unlist()
#vars_shrink_limmar_hom = squeezeVar(sample_vars_hom, counts, robust = TRUE)$var.post
#vars_shrink_limma = squeezeVar(sample_vars, counts-1)$var.post
#vars_shrink_vashr = (vash(sqrt(sample_vars), df = mean(counts)-1)$sd.post)^2

## rrmse
#rrmse_limmar = sqrt(mean((known_vars - vars_shrink_limmar)^2))/sqrt(mean((known_vars - sample_vars)^2))
#rrmse_limma = sqrt(mean((known_vars - vars_shrink_limma)^2))/sqrt(mean((known_vars - sample_vars)^2))
#rrmse_vashr = sqrt(mean((known_vars - vars_shrink_vashr)^2))/sqrt(mean((known_vars - sample_vars)^2))

############### Homoscedastics Scans #############################

# HOMOSCEDASTIC LM
pvalues_lm <- vector("double", length(G))
for (i in 1:(length(G))){
  fit <- lm(Y_hom ~ G[,i] %>% unlist()) %>% summary()
  pvalues_lm[i] <- fit$coefficients[2,4]
}

# HOMOSCEDASTIC EMMA
emm <- emma.ML.LRT(Y_hom, t(G), K)
pvalues_emm <- emm$ps

# HOMOSCEDASTIC KNOWN WISAM
pvalues_wisam_hom <- 
  wisam(G = G, y = YL_hom, strains = t(Z), 
        K = K, weights = "counts")$pvalue

# HOMOSCEDASTIC WITH SAMPLE VARS WEIGHTS 
pvalues_wisam_est_hom <-   
  wisam(G = G, y = YL_hom, strains = t(Z), 
                                 K = K, weights = "samplevars")$pvalue

# NEED HOMOSCEDASTIC WITH SHRUNKEN VARIANCE 
pvalues_wisam_shr_hom <- wisam(G = G, y = YL_hom, strains = t(Z), 
                                 K = K, weights = "limma")$pvalue

# HOMOSCEDASTIC EMMA INDIVIDUAL
emm_hom_ml_z <- emma.ML.LRT(YL_hom, t(G), K, Z)
pvalues_emm_hom_ml_z <- emm_hom_ml_z$ps

print("homosced scans done")

######################## Heteroscedastic Scans ###############################

# HETEROSCEDASTIC LM
pvalues_lm_het <- vector("double", length(G))
for (i in 1:(length(G))){
  fit <- lm(Y ~ G[,i] %>% unlist()) %>% summary()
  pvalues_lm_het[i] <- fit$coefficients[2,4]
}

# HETEROSCEDASTIC EMMA 
emm_het_ml <- emma.ML.LRT(Y, t(G), K)
pvalues_emm_het_ml <- emm_het_ml$ps

# HETEROSCEDASTIC WITH SAMPLE VAR WEIGHTS 
# genome scan with weights
# estimated weights are sample size / std deviation of obs. phenotype per strain
# weights_est = apply(obs_phenotypes_het_tot, 2, function(x) 10/var(x))
pvalues_wisam_est <- wisam(G = G, y = YL, strains = t(Z), 
                             K = K, weights = "samplevars")$pvalue

# HETEROSCEDASTIC WITH VARIANCE SHRINKAGE WEIGHTS 
## variance shrinkage estimation
pvalues_wisam_limmar <- wisam(G = G, y = YL, strains = t(Z), 
                                K = K, weights = "limma")$pvalue

# HETEROSCEDASTIC WITH KNOWN WEIGHTS 
# genome scan with weights
# input known noise/counts instead of sample (i.e. variance as below)
known_weights = counts/variance
known_weights = known_weights/sum(known_weights)*sum(counts)
pvalues_wisam_known <-  wisam(G = G, y = YL, strains = t(Z),
                               K = K, weights = "user", user_weights = known_weights)$pvalue

# HETEROSCEDASTIC EMMA INDIVIDUAL 
emm_het_ml_z <- emma.ML.LRT(YL, t(G), K, Z)
pvalues_emm_het_ml_z <- emm_het_ml_z$ps


## create output tibble 

sim_results<- tibble(d = c(rep(0, SNP-1), 1, rep(0, ncol(G)-SNP)), 
                     # only a 1 where there is the effected SNP 
                              lm = pvalues_lm,
                              emma = pvalues_emm %>% as.vector(),
                              wisam_hom = pvalues_wisam_hom,
                              wisam_est_hom = pvalues_wisam_est_hom,
                              wisam_shr_hom = pvalues_wisam_shr_hom,
                              lm_het = pvalues_lm_het,
                              emma_het_ml = pvalues_emm_het_ml %>% as.vector(),
                              wisam_est = pvalues_wisam_est,
                              wisam_known = pvalues_wisam_known,
                              wisam_limmar = pvalues_wisam_limmar,
			      emma_z_het = pvalues_emm_het_ml_z,
			      emma_z_hom = pvalues_emm_hom_ml_z
)

### roc plot data -- only save the False Positive Rate at the SNP with the true effect 
data = sim_results
output = data.frame(matrix(ncol = 13, nrow = 1))
names(output) = colnames(sim_results)
dist = c(0)#,25,50)
#for(i in 1:3){
    i = 1
    j = which(data$d == 1)
    m =  seq(j-dist[i],j+dist[i])
    pval = apply(data[m,], 2,min)[-1]
    fprs = apply(sweep(data[-m,-1], 2, pval, "-"), 2, function(x) sum(x < 0)/nrow(data))
    output[i,] = c(dist[i], fprs)
#}

### save output 

write_csv(output, output_file)

print("output made")

