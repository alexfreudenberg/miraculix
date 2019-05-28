
.genZuchtwert.orig <- function(Y, X=rep(1, n), G=rep(1,n), Ce=rep(1,n), Cw, maxratio=1000,
                         tol = 0.001, do.plot=FALSE, marginal=NA, trace=TRUE) {  
  cat("entering ", date(), "\n")
  n <- length(Y)
  k <- nrow(Cw)
  stopifnot((is.vector(Ce) && length(Ce) == n) || all(dim(Ce) == n),
            nrow(X) == n,
            (is.vector(G) && length(G) == n) || (nrow(G) == n && ncol(G) == k),
            all(dim(Cw) == k))
  lower <- -(upper <- abs(log(maxratio)))
  if (is.na(marginal)) marginal <- k >= n / 10 

  invCeX <-if (is.vector(Ce))  X / Ce else solvex(Ce, X) # C_e^{-1} X
 # if (is.vector(Ce)) stopifnot( all(X / Ce == solvex(diag(Ce), X)))
  
  XCXXC <- solvex(t(X) %*% invCeX, t(invCeX))
  H <- X %*% XCXXC
  Hstar <- diag(n) - H
  invCstar <- t(Hstar) %*% (if (is.vector(Ce))  Hstar / Ce else solvex(Ce, Hstar))
#  if (is.vector(Ce)) stopifnot(all(Hstar / Ce == solvex(diag(Ce), Hstar)))
  invCw <- solvex(Cw)
  invCstarG <- if (is.vector(G)) invCstar * rep(G, each=n)  else invCstar %*% G
#  if (is.vector(G)) stopifnot(all(invCstar * rep(G, each=n) == invCstar %*% diag(G)))
  GCG <- if (is.vector(G)) G * invCstarG else t(G) %*% invCstarG
#  if (is.vector(G)) stopifnot(all(G * invCstarG == t(diag(G)) %*% invCstarG))
  GCY <- t(invCstarG) %*% Y
  cat("precalculations done ", date(), "\n")
  zaehler <- 0
  ENV <- environment()

  if (marginal) {
    GCwG <- if (is.vector(G)) G * Cw * rep(G, each=n) else G %*% Cw %*% t(G)
  # if (is.vector(G)) stopifnot(all(G * Cw * rep(G, each=n)==diag(G) %*% Cw %*% t(diag(G))))
   ll1 <- function(xi) {
      xi <- exp(xi)
      if (is.vector(Ce)) {
        M <- GCwG
        diag(M) <- diag(M) + Ce * xi
      } else M <- Ce * xi + GCwG
      invMX <- solvex(M, X) # C_e^{-1} X
      beta <- solvex(t(X) %*% invMX, t(invMX) %*% Y)
      Ybeta <- Y - X %*% beta
      vw <- t(Ybeta) %*% solvex(M, Ybeta) / n
      ergb <- n * log(vw) + determinant(M)$modulus
      
      if (trace) {assign("zaehler", zaehler + 1, envir=ENV); cat(zaehler, ": xi=", xi, "  fctn=", formatC(ergb, digits=10), "\n", sep="")}
    return(ergb)
    }
   ll2 <- function(xi) {
 #    print(xi)
 #    str(Ce)
     xi <- exp(xi)
     if (is.vector(Ce)) {
       M <- GCwG
       diag(M) <- diag(M) + Ce * xi        
     } else M <- Ce * xi + GCwG
     invMXY <- solvex(M, cbind(Y, X)) # C_e^{-1} X
     invMX <- invMXY[, -1]
     beta <- solvex(t(X) %*% invMX, t(invMX) %*% Y)
     vw <- t(Y - X %*% beta) %*% (invMXY[, 1] - invMX %*% beta) / n
     ergb <- n * log(vw) + determinant(M)$modulus      
     # stopifnot(all.equal(ergb, ll1(log(xi))))
       if (trace) {assign("zaehler", zaehler + 1, envir=ENV); cat(zaehler, ": xi=", xi, "  fctn=", formatC(ergb, digits=10), "\n", sep="")}
   return(ergb)
    }
    ll <- ll2 ## ll1 is 10% slower
  } else {
    ll <- function(xi) {
      xi <- exp(xi)    
      wxi <- solvex(GCG + xi * invCw, GCY)    
      YGwxi <- Y - (if (is.vector(G)) G * wxi else G %*% wxi)
      Vxi <- t(YGwxi) %*% invCstar %*% YGwxi
      Wxi <- t(wxi) %*% invCw %*% wxi
      ergb <- - k * log(xi) + (k + n) * log(Vxi + xi * Wxi)
     if (trace) {assign("zaehler", zaehler + 1, envir=ENV); cat(zaehler, ": xi=", xi, "  fctn=", formatC(ergb, digits=10), "\n", sep="")}
    return(ergb)
    }
  } 
      
  if (do.plot) {
    m <- ll(0)
#    print(m)
    plot(0, m, xlim=c(lower, upper), ylim=m + c(-1,1) * 20, pch=20,
         main="log-likelihood")
    for (x in seq(0.001, upper, len=100)) {
      min <- ll(-x)
      max <- ll(x)
 #     print(c(-x,x,min,max, is.vector(G), is.vector(Ce)))
      points(c(-x, x),  c(min, max), pch=20)
    }
  }

  res <- optimize(ll, c(lower, upper), tol=tol)
  xi <- exp(res$minimum)

#  print(xi, digits=10)
  
  cat("optimization done ", date(), "\n")
  wxi <- solvex(GCG + xi * invCw, GCY)
  YGwxi <- Y - (if (is.vector(G)) G * wxi else G %*% wxi)
#   if (is.vector(G)) stopifnot(all(G * wxi == diag( G) %*% wxi))

  if (marginal) {
    if (is.vector(Ce)) {
      M <- GCwG
      diag(M) <- diag(M) + Ce * xi
    } else M <- Ce * xi + GCwG    
    invMX <- solvex(M, X) # C_e^{-1} X
    beta0 <- solvex(t(X) %*% invMX, t(invMX) %*% Y)
    Ybeta <- Y - X %*% beta0
    vw <- t(Ybeta) %*% solvex(M, Ybeta) / n
  } else {
    Vxi <- t(YGwxi) %*% invCstar %*% YGwxi
    Wxi <- t(wxi)  %*% invCw %*% wxi
    vw <- (Vxi / xi + Wxi) / (k + n)
  }
  
  beta <- XCXXC %*% YGwxi
  cat("done ", date(), "\n")
  return(list(beta=beta, w=wxi, vw=vw, ve=xi * vw, eps=YGwxi - X %*% beta,
              marginal = marginal, beta0=if (marginal) beta0))
}

