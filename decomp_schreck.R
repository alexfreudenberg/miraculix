## This file is part of the preprint 
#
# Empirical decomposition of the explained variation in the variance components form of the mixed model
# Nicholas Schreck 2019
#
## We include it for documentary purposes 


#Packages
library(sommer)
library(BGLR)
################################################################################
### Data
data(mice)
X <- scale(mice.pheno$Obesity.BodyLength, center=TRUE, scale=FALSE)
n <- length(X)
Z <- mice.X
p <- dim(Z)[2]
y <- mice.pheno$Obesity.BMI
rm(mice.A, mice.pheno, mice.X)
decomp <- Var_decomp_lmm(y=y, X=X, Z=Z) # Result for the "Real Data Application"


################################################################################
### Functions
Var_decomp_lmm <- function(y, X, Z){
  ##############################################################################
  # Input: Phenotypic vector, Matrix of fixed effects, SNP genotype matrix
  # Aim  : Calculate variance decompositon in the LMM in the var. comp. form
  # Output : Components of the variance decomposition
  ##############################################################################
  ### Data handling
  n <- length(X[,1])
  k <- dim(X)[2]
  y.scaled <- y / sd(y)
  Zc <- scale(Z, center=T, scale=FALSE)
  G <- tcrossprod(Zc)
  rm(Zc)
  Xc <- scale(X, center=TRUE, scale=FALSE)
  Sx <- crossprod(Xc) / (n - 1)
  ### Fit LMM in var. comp. form (equivalent parameterization)
  rownames(G) <- as.character(1:n)
  colnames(G) <- as.character(1:n)
  cat("Start MME")
  library(sommer)
  lmm <- mmer(Y~1+Xc, random=~vs(factor(1:n), Gu=G), method = "AI",
              data=data.frame(Y=y.scaled, Xc=Xc), verbose=T, tolpar=1e-6,
              date.warning=FALSE)
  ### Calculate variance decomposition
  var.g <- as.numeric(lmm$sigma$`u:1:n`) # variance component g
  var.e <- as.numeric(lmm$sigma$units) # variance component eps
  gv_gcta <- var.g * sum( diag(G) ) / (n - 1)
  b <- as.vector(lmm$Beta$Estimate) # fixed effect estimate
  g <- as.vector(lmm$U$`u:1:n`$Y) # empirical best predictor of g=Zu
  cov.b <- as.matrix(lmm$VarBeta) # covariance matrix of BLUE
  R2x <- t(b[-1]) %*% Sx %*% b[-1] - 
    sum( diag( Sx %*% cov.b[2:(k+1), 2:(k+1)]) )
  R2xz <-  2 * t(b[-1]) %*% crossprod(Xc, g) / (n - 1)
  #  rm(Sx, Xc, b, cov.b)
  cov.g <- as.matrix(lmm$VarU$`u:1:n`$Y) 
  #  rm(lmm)
  R2z <- gv_gcta + sum(g ^ 2) / (n - 1) - sum( diag(cov.g) ) / (n - 1)
  # rm(cov.g, G)
  summe <- R2x + R2z + R2xz + var.e
  ### Return
  return(list(Rx=R2x, Rz=R2z, Rxz=R2xz, vare=var.e, varg=var.g, gcta=gv_gcta,
              total=summe)  )
}
################################################################################