

.genZuchtwert <- function(Y, X=rep(1, n), G=rep(1,n), Ce=rep(1,n), Cw, maxratio=1000,
                         tol = 0.001, do.plot=FALSE) {
  # crossproduct sollte schneller sein. Ist es aber nicht!
 cat("entering ", date(), "\n")
  n <- length(Y)
  k <- nrow(Cw)
  stopifnot(all(dim(Ce) == n), nrow(X) == n, nrow(G) == n, ncol(G) == k,
            all(dim(Cw) == k))
  invCeX <- solvex(Ce, X) # C_e^{-1} X
  XCXXC <- solvex(crossprod(X, invCeX), t(invCeX))
  H <- X %*% XCXXC
  Hstar <- diag(n) - H
  cat("H calculated ", date(), "\n")
  invCstar <- crossprod(Hstar, solvex(Ce, Hstar))
  cat("invCwstar ", date(), "\n")
  invCw <- solvex(Cw)
  cat("invCw ", date(), "\n")
  invCstarG <- invCstar %*% G
  cat("invCstarG ", date())
  GCG <- crossprod(G, invCstarG)
  cat("GCG ", date(), "\n")
  GCY <- crossprod(invCstarG, Y)
  cat("GCY ", date(), "\n")

  ll <- function(xi) {
    #print(xi)
    xi <- exp(xi)    
    wxi <- solvex(GCG + xi * invCw, GCY)
    YGwxi <- Y - G %*% wxi
    Vxi <- crossprod(YGwxi, invCstar %*% YGwxi)
    Wxi <- crossprod(wxi, invCw %*% wxi)
    return(- k * log(xi) + (k + n) * log(Vxi + xi * Wxi))
  }

  lower <- -(upper <- abs(log(maxratio)))
#  print(c(lower, upper))
 
  if (do.plot) {
    m <- ll(0)
 #   print(m)
    plot(Inf, Inf, xlim=c(lower, upper), ylim=m + c(-1,1) * 2000)
    for (x in seq(0.001, upper, len=100)) {
      min <- ll(-x)
      max <- ll(x)
  #    print(c(-x,x,min,max, is.vector(G), is.vector(Ce)))
      points(c(-x, x),  c(min, max), pch=20)
    }
  }


 res <- optimize(ll, c(lower, upper), tol=tol)
  xi <- exp(res$minimum)
  
  wxi <- solvex(GCG + xi * invCw, GCY)
  YGwxi <- Y - G %*% wxi
  Vxi <- t(YGwxi) %*% invCstar %*% YGwxi
  Wxi <- t(wxi)  %*% invCw %*% wxi
  vw <- (Vxi / xi + Wxi) / (k + n)
  
  beta <- XCXXC %*% YGwxi
  cat("done ", date(), "\n")
  return(list(beta=beta, w=wxi, vw=vw, ve=xi * vw, eps=YGwxi))
}
