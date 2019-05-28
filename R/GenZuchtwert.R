
genZuchtwert <- function(Y, X=rep(1, n), G=rep(1,n), Ce=rep(1,n), Cw, maxratio=1000,
                         tol = 0.001, do.plot=FALSE, marginal=NA, REML=FALSE, trace=TRUE) {
  cat("entering ", date(), "\n")
  n <- length(Y)
  k <- nrow(Cw)
  stopifnot((is.vector(Ce) && length(Ce) == n) || all(dim(Ce) == n),
            nrow(X) == n,
            (is.vector(G) && length(G) == n) || (nrow(G) == n && ncol(G) == k),
            all(dim(Cw) == k))
  if (is.na(marginal)) marginal <- k >= n / 10
  if (REML) .genZuchtwertREML(Y=Y, X=X, G=G, Ce=Ce, Cw=Cw, maxratio=maxratio,
                             tol=tol, do.plot=do.plot, trace=trace)
  else if (marginal) .genZuchtwertMarginal(Y=Y, X=X, G=G, Ce=Ce, Cw=Cw, maxratio=maxratio,
                                          tol=tol, do.plot=do.plot, trace=trace)
  else .genZuchtwertFull(Y=Y, X=X, G=G, Ce=Ce, Cw=Cw, maxratio=maxratio,
                        tol=tol, do.plot=do.plot, trace=trace)
}





.genZuchtwertMarginal <- function(Y, X=rep(1, n), G=rep(1,n), Ce=rep(1,n), Cw,
                                  maxratio=1000, tol = 0.001, do.plot=FALSE, trace=TRUE) {
  n <- length(Y)
  k <- nrow(Cw)
  lower <- -(upper <- abs(log(maxratio)))
  CwtG <- if (is.vector(G)) Cw * rep(G, each=n) else  Cw %*% t(G) 
  GCwG <- if (is.vector(G)) G * CwtG else G %*% CwtG
  diagGCwG <- diag(GCwG) 
  ##  GCwG <- if (is.vector(G)) G * Cw * rep(G, each=n) else G %*% Cw %*% t(G)
  cat("precalculations done ", date(), "\n")
  zaehler <- 0
  ENV <- environment()
  beta <- MYXb <- vw <- xi <- NULL

  ll <- function(xi) {
    assign("xi", xi <- exp(xi), envir=ENV)
    if (is.vector(Ce)) {
      M <- GCwG
      diag(M) <- diagGCwG + Ce * xi        
    } else M <- Ce * xi + GCwG  
    invdet <- solvex(M, cbind(Y, X), logdeterminant = TRUE) # C_e^{-1} X
    invMX <- invdet$inv[, -1]
    assign("beta", solvex(t(X) %*% invMX, t(invMX) %*% Y), envir=ENV)
    assign("MYXb", invdet$inv[, 1] - invMX %*% beta, envir=ENV) # M^{-1}(Y - Xb)
    assign("vw", ( t(Y - X %*% beta) %*% MYXb ) / n, envir=ENV)
    ergb <- n * log(vw) + invdet$logdet
    if (trace) {
      assign("zaehler", zaehler + 1, envir=ENV);
      cat(zaehler, ". xi=", xi, "  fctn=", formatC(ergb, digits=10), "\n",
          sep="")
    }
    return(ergb)
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
  cat("optimization done ", date(), "\n")
  
  if (abs(exp(res$minimum) - xi) > 1e-8) ll(res$minimum)
  wxi <- CwtG %*% MYXb
  YGwxi <- Y - (if (is.vector(G)) G * wxi else G %*% wxi)
  ##   if (is.vector(G)) stopifnot(all(G * wxi == diag( G) %*% wxi)) 

  cat("done ", date(), "\n")
  return(list(beta=beta, w=wxi, vw=vw, ve=xi * vw, eps=YGwxi - X %*% beta, marginal = TRUE))
}



.genZuchtwertREML <- function(Y, X=rep(1, n), G=rep(1,n), Ce=rep(1,n), Cw,
                             maxratio=1000,
                             tol = 0.001, do.plot=FALSE, trace=TRUE) {
  n <- length(Y)
  k <- nrow(Cw)
  lower <- -(upper <- abs(log(maxratio)))
  CwtG <- if (is.vector(G)) Cw * rep(G, each=n) else  Cw %*% t(G) 
  GCwG <- if (is.vector(G)) G * CwtG else G %*% CwtG
  ##  GCwG <- if (is.vector(G)) G * Cw * rep(G, each=n) else G %*% Cw %*% t(G)
  diagGCwG <- diag(GCwG)
  XXtX <- solvex(t(X) %*% X, t(X))
  H <- (diag(n) - X %*% XXtX) [-n , ]
  HY <- H %*% Y
  HGCwGH <- H %*% GCwG %*% t(H)
  
  HCeH <- H %*% (if (is.vector(Ce))  Ce * t(H) else Ce %*% t(H))

  
  cat("precalculations done ", date(), "\n")
  zaehler <- 0
  ENV <- environment()
  beta <- MYXb <- vw <- xi <- NULL

  ll <- function(xi) {
    assign("xi", xi <- exp(xi), envir=ENV)
    M <- HCeH * xi + HGCwGH
    invdet <- solvex(M, HY, logdeterminant=TRUE)
    assign("vw", ( t(HY) %*% invdet$inv) / (n-1), envir=ENV)
    ergb <- (n - 1) * log(vw) + invdet$logdet   
    if (trace) {
      assign("zaehler", zaehler + 1, envir=ENV);
      cat(zaehler, ". xi=", xi, "  fctn=", formatC(ergb, digits=10), "\n", sep="")
    }
    return(ergb)
  }

  llorig <- function(xi) {
    assign("xi", xi <- exp(xi), envir=ENV)
    if (is.vector(Ce)) {
      M <- GCwG
      diag(M) <- diagGCwG + Ce * xi        
    } else M <- Ce * xi + GCwG
    invdet <- solvex(M, cbind(Y, X), logdeterminant = TRUE) # C_e^{-1} X
    invMX <- invdet$inv[, -1]
    assign("beta", solvex(t(X) %*% invMX, t(invMX) %*% Y), envir=ENV)
    assign("MYXb", invdet$inv[, 1] - invMX %*% beta, envir=ENV) # M^{-1}(Y - Xb)
    assign("vw", ( t(Y - X %*% beta) %*% MYXb ) / n, envir=ENV)
    ergb <- n * log(vw) + invdet$logdet
    if (trace) {
      assign("zaehler", zaehler + 1, envir=ENV);
      cat(zaehler, ". xi=", xi, "  fctn=", formatC(ergb, digits=10), "\n", sep="")
    }
    return(ergb)
  }
 
  
  if (do.plot) {
    m <- ll(0)
#    print(m)
    plot(0, m, xlim=c(lower, upper), ylim=m + c(-1,1) * 20, pch=20,
         main="log-likelihood")
    for (x in seq(0.001, upper, len=100)) {
      min <- ll(-x)
      max <- ll(x)
#      print(c(-x,x,min,max, is.vector(G), is.vector(Ce)))
      points(c(-x, x),  c(min, max), pch=20)
    }
  }

  boundary <- c(lower, upper)
  res <- optimize(ll, boundary, tol=tol)
  cat("optimization done ", date(), "\n")
  
  llorig(res$minimum)
  wxi <- CwtG %*% MYXb
  YGwxi <- Y - (if (is.vector(G)) G * wxi else G %*% wxi)
  ##   if (is.vector(G)) stopifnot(all(G * wxi == diag( G) %*% wxi)) 

  cat("done ", date(), "\n")
  return(list(beta=beta, w=wxi, vw=vw, ve=xi * vw, eps=YGwxi - X %*% beta,
              boundary=boundary,
              marginal = TRUE))
}



.genZuchtwertFull <- function(Y, X=rep(1, n), G=rep(1,n), Ce=rep(1,n), Cw, maxratio=1000,
                         tol = 0.001, do.plot=FALSE, trace=TRUE) {
  marginal <- FALSE
  n <- length(Y)
  k <- nrow(Cw)
   lower <- -(upper <- abs(log(maxratio)))
  if (is.na(marginal)) marginal <- k >= n / 10 

  invCeX <- if (is.vector(Ce))  X / Ce else solvex(Ce, X) # C_e^{-1} X
 # if (is.vector(Ce)) stopifnot( all(X / Ce == solvex(diag(Ce), X)))
  
  XCXXC <- solvex(t(X) %*% invCeX, t(invCeX))
  H <- X %*% XCXXC
  Hstar <- diag(n) - H
  invCstar <- t(Hstar) %*% (if (is.vector(Ce))  Hstar / Ce
                            else solvex(Ce, Hstar))
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

  ll <- function(xi) {
    xi <- exp(xi)    
    wxi <- solvex(GCG + xi * invCw, GCY)    
    YGwxi <- Y - (if (is.vector(G)) G * wxi else G %*% wxi)
    Vxi <- t(YGwxi) %*% invCstar %*% YGwxi
    Wxi <- t(wxi) %*% invCw %*% wxi
    ergb <- - k * log(xi) + (k + n) * log(Vxi + xi * Wxi)
    if (trace) {assign("zaehler", zaehler + 1, envir=ENV); cat(zaehler, ". xi=", xi, "  fctn=", formatC(ergb, digits=10), "\n", sep="")}
    return(ergb)
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

  Vxi <- t(YGwxi) %*% invCstar %*% YGwxi
  Wxi <- t(wxi)  %*% invCw %*% wxi
  vw <- (Vxi / xi + Wxi) / (k + n)
  
  beta <- XCXXC %*% YGwxi
  cat("done ", date(), "\n")
  return(list(beta=beta, w=wxi, vw=vw, ve=xi * vw, eps=YGwxi - X %*% beta,
              marginal = marginal))
}

