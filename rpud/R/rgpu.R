# Author Alexander Freudenberg
# Copyright 2021 Alexander Freudenberg

cu_matmul <- function(A,B,operation=0L){
	.Call(C_cu_matmul,A,B,as.integer(operation))
}	
	
access_pointer <- function(ptr){
  .Call(C_access_pointer, ptr)
}

solve_large <- function(M,b){
  # Calculates the solution to large equation systems
  if(nrow(M)!=ncol(M))error("Wrong dims")
  n<- nrow(M)
  n1 <- as.integer(n/2)
  n2 <- n-n1

  # Creates block matrices
  A<- M[1:n1, 1:n1]
  B<- M[1:n1, (n1+1):n]
  C<- M[(n1+1):n, 1:n1]
  D<- M[(n1+1):n, (n1+1):n]
  b1<- as.matrix(b[1:n1,])
  b2<- as.matrix(b[(n1+1):n,])

  # Solves sub-equation systems
  print(system.time(AinvBb <- solvex(A,cbind(B,b1))))
  Y <- AinvBb[,1:n2]
  y_2 <- AinvBb[,(n2+1): (n2+length(b1[1,])) ]
  print(system.time(S <- D -cu_matmul( t(C), Y, 0L)))
  print(system.time(y_4<- solvex(S, C %*% y_2 - b2)))

  # Calculates final result
  x<- rbind(y_2 + Y %*% y_4, -y_4)
  return(x)
}

cu_decomp <- function(y,X,Z){
  # Calculates the empirical decomposition of the phenotypic variance for the mixed effects model
  n<-length(y)
  stopifnot(dim(X)[1]==length(y)&& dim(Z)[1]==length(y))
  k<-dim(X)[2]
  y.scaled <- y/sd(y)

  # Calculate GRM and sample covariance matrix
  G<- miraculix::relationshipMatrix(t(Z),centered=T,normalized=F)
  Sigma_X <- crossprod(scale(X,center=TRUE,scale = FALSE))/(n-1)

  # Determine REML estimates
  lmm <- solveMME(X,G,y.scaled,start_params = c(1,1))
  var.g<- lmm$var.g
  var.epsilon <- lmm$var.e
  gv_gcta <- var.g * sum(diag(G))/(n-1)
  beta <- lmm$Betahat[-1]
  g <- lmm$g
  cov_beta <- lmm$CovBeta[-1,-1]
  cov_g <- lmm$Covg 
  
  # Determine REML estimates of the empirical variances
  R2x <-  t(beta)%*% crossprod(Sigma_X[2,2],beta) - 
    sum(diag(cu_matmul(Sigma_X[2,2], cov_beta)))
  R2xz <- 2 * (beta )%*% crossprod(X[,2], g-mean(g))/(n-1)
  R2z <- gv_gcta + sum(g^2)/(n-1)-sum(diag(cov_g))/(n-1)

  return(list(R2x=R2x, R2z=R2z, R2xz=R2xz, var.e= var.epsilon, var.g=var.g, gcta=gv_gcta, total=R2xz+R2z+R2x + var.epsilon))
}

solveMME <- function(X,G,y, start_params=c(1,1),maxiter=20,toler=1e-3){
  # Numerically determines the REML estimates for variance components
  n<-length(y)
  m <- n - dim(X)[2]

  # Calculate the logdet
  logdet_XtX <- solvex(cu_matmul(X,X),logdeterminant = TRUE)$logdet
  params <- start_params
  print(params)

  # Initialize parameters
  H<- build_H(G,n,params)
  P<- 1/params[1] * build_P(H,X,n)
  P <- as.numeric(m/(t(y)%*%P%*%y))*P
  PG <-  cu_matmul(P,G)
  GPy <- t(PG) %*% y
  Py <- P%*%y

  # Determine gradient and likelihood
  grad <- eval_grad(P,PG,Py,GPy,params)
  llik_new <- eval_likelihood(Py,y,1/params[1]*H,X, logdet_XtX)
  llik <- Inf
  print(grad)
  iter <- 1
  
  # Calculate first updates with Fisher scoring algorithm
  while(abs(llik-llik_new)>1&& iter<maxiter){
 # for(iter in 1:3){
    if(sum(grad^2)<toler)break;
    # Calculate Fisher information matrix
      M <- matrix(ncol=2,nrow=2)
    M[1,1] <- 0.5*sum(diag(cu_matmul(P,P)))
    M[1,2] <- M[2,1] <- 0.5*sum(diag(cu_matmul(t(PG),P)))
    M[2,2] <-  0.5*sum(diag(cu_matmul(t(PG),PG)))
    if(det(M)<1e-2)break;
    # Update parameters
    params <- params +0.7*solvex(M,grad)
    params[params<0] <- 1e-3
    print(params)

    #Update matrices
    H<- build_H(G,n,params)
    P<- 1/params[1] * build_P(H,X,n)
    rss <- as.numeric( m/(t(y)%*%P%*%y) )
    print(rss)
    P <-rss *P
    params <- params/rss
    PG <-  cu_matmul(P,G)
    GPy <- t(PG) %*% y
    Py <- P%*%y
    
    #Caclulate gradient, likelihood
    grad <- eval_grad(P,PG,Py,GPy,params)
    print(grad)
    llik <- llik_new
    llik_new <- eval_likelihood(Py,y,1/params[1]*H,X, logdet_XtX)
    #if(abs(llik_new-llik)<toler)break;
    cat(paste("Iter:", iter, ", Likelihood:", llik_new,", Grad:",sum(grad^2), ", Params:",paste(c(params[1],params[2]),collapse = ",")))
    iter <- iter+1
  }
  # Switch to Average Information
 cat("\nSwitch\n")
  while(sum(grad^2)>toler &&iter <maxiter){
    #Calculate Average Information matrix
    M <- calc_AI(P,Py,GPy,y,params)

  # Stretch if singular
    if(abs(det(M))<1e-1)M <- M+1e-1*diag(2)
    params <- params +0.9*solvex(M,grad)
    params[params<0] <- 1e-3

    #Update parameters
    H <- build_H(G,n,params)
    P <- 1/params[1]*build_P(H,X,n)
    rss <- as.numeric( m/(t(y)%*%P%*%y) )
    print(rss)
    P <-rss*P
    params <- params/rss
    PG <-  cu_matmul(P,G)
    GPy <- t(PG) %*% y
    Py <- P%*%y

    #Calculate gradient, likelihood
    grad <- eval_grad(P,PG,Py,GPy,params[1])
    llik_new <- eval_likelihood(Py,y,1/params[1]*H,X, logdet_XtX)
    llik <- llik_new
    cat(paste("Iter:", iter, ", Likelihood:", llik,", Grad:",sum(grad^2), ", Params:",paste(c(params[1],params[2]),collapse = ",")))
    iter<-iter+1
  }
  #Calculate output data
  Sigma_xi <- params[1]/params[2]*H
  SigmaxiX <- solvex(Sigma_xi,X)
  betahat <- solvex(cu_matmul(X,SigmaxiX)) %*% (t(SigmaxiX)%*%y)
  g <-  G %*% solvex(Sigma_xi, (y-X %*%betahat))
  XHiXt_i <- solve(cu_matmul(X,solvex(H,X)))
  CovBeta <- XHiXt_i*params[1]
  HiG <- solvex(H,G)
  X_XHiXt_i_Xt <- crossprod(t(X),crossprod(t(XHiXt_i),t(X)))
  Covg <- params[2]^2/params[1] * (cu_matmul(t(HiG),G) - cu_matmul(cu_matmul(HiG, X_XHiXt_i_Xt),HiG))
  return(list(var.e=params[1],var.g = params[2], g=g, Betahat = betahat, CovBeta = CovBeta, Covg = Covg, llik = llik, converged= sum(grad^2)<toler ))
}

build_H <- function(G,n,params) diag(n) + params[2]/params[1] * G
build_P <- function(H,X,n){
  # #Buffer matrix
  A <- solvex(H, cbind(diag(n),X))
  VinvX <- A[,(n+1):(n+dim(X)[2])]
  #VińvX <- solvex(H,X)
  XVinvX_inv <- solvex( cu_matmul(X, VinvX), t(VinvX) )
   P <- A[,1:n] - crossprod(t(VinvX), XVinvX_inv)
  #P <- solvex(H,diag(n)-crossprod(t(X),XVinvX_inv) )
  return(P)
}
eval_grad <- function(P,PG,Py,GPy,params){
  grad_e <- -1/2* sum(diag(P))  + 1/2 * t(Py)%*%Py 
  grad_k <- -1/2* sum(diag( PG )) + 1/2* ( t(Py) )%*% (GPy)
  return(c(grad_e,grad_k))
}
eval_likelihood <- function(Py,y,V,X, logdet_XtX){
  V_list <- solvex(V,X,logdeterminant = T)
  logdet_sum <- -logdet_XtX + V_list$logdet + solvex(cu_matmul(X,V_list$inv),logdeterminant = TRUE)$logdet
  return(1/2*logdet_sum - 1/2 * t(y)%*%Py)
}
calc_AI <- function(P,Py,GPy,y,params){
  M<- matrix(ncol= length(params),nrow=length(params))
  M[1,1] <- 1/(2)* (t(Py)%*%P%*%P%*%P%*% Py)
  M[1,2] <- M[2,1] <- 1/2* t(Py)%*%P%*%P%*%P%*% GPy 
  M[2,2] <- 1/2 * t(GPy)%*%P%*%P%*%P%*% GPy
  return(M)
}
