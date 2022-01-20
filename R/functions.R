# this function is used to find the h2 that optimizes the log likelihood
# h2 is the value that is changed to optimize the log likelihood
# X is matrix of covariates + SNP
# y is phenotype vector
# K is genetic relationship matrix
# Dhalf is a function of the weight matrix (diagonal matrix of square root of weights)
# p is the number of strains (make this consistent)
# laml is the eigenvalues of L
# Ul is the eigenvectors of L
# I is a p x p identity
likeli.brent.het <- function(h2, X, y, K, Dhalf, p, laml, Ul, I){
  Mhet = diag(1/sqrt(h2*laml + (1-h2)*diag(I))) %*% t(Ul) %*% Dhalf
  Mx = Mhet %*% X
  My = Mhet %*% y
  fit = stats::lm.fit(Mx, c(My))
  B = fit$coefficients
  s2 = sum(fit$residuals^2)/p

  # log determinant V
  LDV = sum(log(diag(Dhalf))) + sum(log(h2*laml + (1-h2)))

  # log likelihood
  log = -.5*LDV + (-p/2)*log(s2) + (-p/2)*log(2*pi) + (-.5*(p))
  log
}

likeli.brent.het.estimates <- function(h2, X, y, K, Dhalf, p, laml, Ul, I){
  Mhet = diag(1/sqrt(h2*laml + (1-h2)*diag(I))) %*% t(Ul) %*% Dhalf
  Mx = Mhet %*% X
  My = Mhet %*% y
  fit = stats::lm.fit(Mx, c(My))
  B = fit$coefficients
  s2 = sum(fit$residuals^2)/p

  list(B = B, s2 = s2)
}


# this function runs a genome scan
# G is genotype matrix
# y is phenotype vector
# X is covariates (or simply, the intercept)
# K is genetic relationship matrix
# weight is a vector of the weights
scan.strain.means <- function(G, y, X, K, weights){
  # p = number of strains
  p = dim(K)[1]
  I = diag(1,p)

  Dhalf = diag(sqrt(weights))

  p_value_het = vector("numeric", length(G)) # p values
  thetaMLEs = vector("numeric", length(G)) # likelihood under alt
  theta0s = vector("numeric", length(G)) # likelihood under null
  s20s = vector("numeric", length(G)) # s2 estimate under null
  h20s = vector("numeric", length(G)) # h2 estimate under null
  betas = vector("numeric", length(G)) # beta estimate under null

  # intercept covariate
  Xint = X
  # L <- t(t(sqrt(w) * K) * sqrt(w)) --- this is robert's way of finding L, faster?
  L = Dhalf %*% K %*% Dhalf
  eL = eigen(L)
  laml = eL$values
  Ul = eL$vectors
  # I = diag(1, p)

  ### RUNNING THE SCAN

  # null scan
  opt0 = optimize(likeli.brent.het, c(0,1), Xint, y, K, Dhalf, p, laml, Ul, I, maximum = TRUE)
  theta0 = opt0$objective
  h20 = opt0$maximum

  #### Extract the parameters
  params = likeli.brent.het.estimates(h20, Xint, y, K, Dhalf, p, laml, Ul, I)
  s20 = params$s2 %>% unname() %>% as.numeric()
  #tau2 = s2*h2_max
  #sig2 = s2-tau2


  # run the single scan for all columns of G
  #
  for (i in 1:length(G)){
    SNP = suppressWarnings(
      strsplit(str_replace_all(G[i], c("NA"="P", "0.5"="5")),
               split ="")[[1]]%>%as.numeric())
    SNP[which(SNP==5)] = 0.5
    # bind intercept with SNP genotypes
    X = cbind(Xint, SNP)
    # check for NAs
    if (any(is.na(X))){
      # take out the NA
      X0 = na.omit(X)
      # find the indices associated with the NA
      index = na.omit(X) %>% na.action()
      # remove those indices from K, y, weights
      K0 <- K[-c(index), -c(index)]
      y0 <- y[-c(index)]
      Dhalf0 <- Dhalf[-c(index), -c(index)]
      # new number of strains
      p0 = dim(K0)[1]
      I0 = diag(1, p0)
      # re-decompose matrices
      L0 = Dhalf0 %*% K0 %*% Dhalf0
      #L0 = (Dhalf0 %*% K0 %*% Dhalf0) %*% t(Z0)
      eL0 = eigen(L0)
      laml0 = eL0$values
      Ul0 =  eL0$vectors
      # new null scan and SNP scan
      opt00 = optimize(likeli.brent.het, c(0,1), X0[,1], y0, K0, Dhalf0, p0, laml0, Ul0, I0,
                       maximum = TRUE)
      theta00 = opt00$objective
      h200 = opt00$maximum

      #### Extract the parameters
      params0 = likeli.brent.het.estimates(h200, X0[,1], y0, K0, Dhalf0, p0, laml0, Ul0, I0)
      beta000 = params0$B %>% unname() %>% as.numeric()
      s200 = params0$s2 %>% unname() %>% as.numeric()
      #tau2 = s2*h2_max
      #sig2 = s2-tau2

      optMLE = optimize(likeli.brent.het, c(0,1), X0, y0, K0, Dhalf0, p0, laml0, Ul0, I0,
                        maximum = TRUE)
      thetaMLE = optMLE$objective
      h2 = optMLE$maximum

      #### Extract the parameters
      paramsMLE = likeli.brent.het.estimates(h2, X0, y0, K0, Dhalf0, p0, laml0, Ul0, I0)
      beta = paramsMLE$B %>% unname() %>% as.numeric()
      beta10 = beta[2]

      tLR = 2*(thetaMLE - theta00)
      thetaMLEs[i] = thetaMLE
      theta0s[i] = theta00
      s20s[i] = s200
      h20s[i] = h200
      betas[i] = beta10
      p_value_het[i] = pchisq(tLR, 1, lower.tail = FALSE)
    } else {
      # SNP scan for no NA
      optMLE = optimize(likeli.brent.het, c(0,1), X, y, K, Dhalf, p, laml, Ul, I,
                        maximum = TRUE)
      thetaMLE = optMLE$objective
      h2 = optMLE$maximum

      #### Extract the parameters
      paramsMLE = likeli.brent.het.estimates(h2, X, y, K, Dhalf, p, laml, Ul, I)
      beta = paramsMLE$B %>% unname() %>% as.numeric()
      beta1 = beta[2]

      tLR = 2*(thetaMLE - theta0)
      thetaMLEs[i] = thetaMLE
      theta0s[i] = theta0
      s20s[i] = s20
      h20s[i] = h20
      betas[i] = beta1
      p_value_het[i] = pchisq(tLR, 1, lower.tail = FALSE)
    }
  }
  # returns p values
  list(ps = p_value_het, theta1 = thetaMLEs, theta0 = theta0s,
       h20 = h20s, s20 = s20s, beta = betas)
}

### used to find unique snps in genome scan
mapping = function(x){
  i = 1
  while(is.na(x[i])|x[i]==0.5){
    i = i+1
  }
  if(x[i]==1){
    x = plyr::mapvalues(x, c(0,1), c(1,0), warn_missing = FALSE)
  }
  x
}

### not written myself: from http://mouse.cs.ucla.edu/emma/install.html
emma.kinship <- function(snps, method="additive", use="all") {
n0 <- sum(snps==0,na.rm=TRUE)
nh <- sum(snps==0.5,na.rm=TRUE)
n1 <- sum(snps==1,na.rm=TRUE)
nNA <- sum(is.na(snps))

stopifnot(n0+nh+n1+nNA == length(snps))

if ( method == "dominant" ) {
  flags <- matrix(as.double(rowMeans(snps,na.rm=TRUE) > 0.5),nrow(snps),ncol(snps))
  snps[!is.na(snps) & (snps == 0.5)] <- flags[!is.na(snps) & (snps == 0.5)]
}
else if ( method == "recessive" ) {
  flags <- matrix(as.double(rowMeans(snps,na.rm=TRUE) < 0.5),nrow(snps),ncol(snps))
  snps[!is.na(snps) & (snps == 0.5)] <- flags[!is.na(snps) & (snps == 0.5)]
}
else if ( ( method == "additive" ) && ( nh > 0 ) ) {
  dsnps <- snps
  rsnps <- snps
  flags <- matrix(as.double(rowMeans(snps,na.rm=TRUE) > 0.5),nrow(snps),ncol(snps))
  dsnps[!is.na(snps) & (snps==0.5)] <- flags[!is.na(snps) & (snps==0.5)]
  flags <- matrix(as.double(rowMeans(snps,na.rm=TRUE) < 0.5),nrow(snps),ncol(snps))
  rsnps[!is.na(snps) & (snps==0.5)] <- flags[!is.na(snps) & (snps==0.5)]
  snps <- rbind(dsnps,rsnps)
}

if ( use == "all" ) {
  mafs <- matrix(rowMeans(snps,na.rm=TRUE),nrow(snps),ncol(snps))
  snps[is.na(snps)] <- mafs[is.na(snps)]
}
else if ( use == "complete.obs" ) {
  snps <- snps[rowSums(is.na(snps))==0,]
}

n <- ncol(snps)
K <- matrix(nrow=n,ncol=n)
diag(K) <- 1

for(i in 2:n) {
  for(j in 1:(i-1)) {
    x <- snps[,i]*snps[,j] + (1-snps[,i])*(1-snps[,j])
    K[i,j] <- sum(x,na.rm=TRUE)/sum(!is.na(x))
    K[j,i] <- K[i,j]
  }
}
return(K)
}




