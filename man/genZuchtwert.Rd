\name{GenZuchtwert}
\alias{genZuchtwert}
\alias{genZuchtwert2}
\title{
 GenZuchtwert 
}
\description{
  We solve the mixed model
  \deqn{ Y = Xb + Gw + \varepsilon}
  where \eqn{b} is a fixed effect and \eqn{w} is a random effect.
  It is assumed that
  \eqn{w} is \eqn{N(0, v_w C_w)} and
  \eqn{\varepsilon} is \eqn{N(0, v_\varepsilon C_\varepsilon)}

  See \href{../doc/genZuchtwert.pdf}{genZuchtwert.pdf} for the
  derivation of the formulae.
}
\usage{
genZuchtwert(Y, X=rep(1, n), G=rep(1,n), Ce=rep(1,n), Cw, maxratio=1000,
             tol = 0.001, do.plot=FALSE, marginal=NA, REML=FALSE,
             trace=TRUE)
}
\arguments{ 
  \item{Y,X}{see above
  }
  \item{G}{vector or matrix. If vector it is interpreted as a
    diagonal matrix
    }
  \item{Ce}{vector or matrix. \eqn{C_\varepsilon} above.
    if vector than it is interpreted as a diagonal matrix.
  }
  \item{Cw}{\eqn{C_w} above
  }
  \item{maxratio}{positive number
    maximal (or minimal) ratio between \eqn{v_w} and \eqn{v_\varepsilon}
    up to (down to) which the algorithm searches for the argmax of
    the likelihood.
  }
  \item{tol}{
    see \link[stats]{optimize}
  }
  \item{do.plot}{
    logical. If \code{TRUE} the likelihood as a function of
    \code{ve}/\code{vw} is plotted before the optimization is performed.
  }
  \item{marginal}{logical or \code{NA}.
    If \code{TRUE} then the overall distribution of \code{Y} is
    considered:
    \eqn{Y \sim N(X b, v_\varepsilon C_\varepsilon + v_w G C_w G^T )}

    If \code{FALSE} then \code{Y} is modelled by a conditional
    distribution given the \code{w} and the distribution \code{w} itself.

    If \code{NA} then \code{marginal=TRUE} iff
    \code{nrow(Cw) > length(Y) / 10}.
  }
  \item{REML}{logical. If \code{TRUE} then REML instead of MLE is
    performed. The argument \code{marginal} given by the user is ignored then (and is
    always considered to be \code{TRUE})
    }
  \item{trace}{logical. If \code{TRUE} then the optimizer is traced.}
}

%\details{}
% http://stackoverflow.com/questions/15397390/parallelized-optimization-in-r

 
\value{
  The function calculates the MLE and returns a list that contains
  \code{beta} (\eqn{\beta}),
  \code{vw} (\eqn{v_w}),
  \code{ve} (\eqn{v_\varepsilon}),
  the actually used value for \code{marginal} 
  and the residuals \code{eps}.
  
  If \code{marginal} is \code{TRUE} a control variable (\code{beta0})
  is returned.
}


\author{Martin Schlather, \email{schlather@math.uni-mannheim.de}
}

\examples{
\dontrun{

n <- 1600 # up to 6000
k <- as.integer(n) # 0.2, 0.1
k <- as.integer(n) # 0.2, 0.1
ve <- .10
Ce <- diag(n)
#
Ce <- rep(1, n); ## if deleted than slower
G <- if (k==n) diag(n) else matrix(nr = n, nc=k, rexp(n * k))
#
G <- rep(1, n) ## if deleted than slower
Cw <- matrix(nr=k, nc=k, rexp(k * k))
Cw <-  Cw \%*\% t(Cw)
vw <- 0.1
w <- t(chol(vw * Cw)) \%*\%  rnorm(k)
mu <- 0.5
eps <-
  if (is.vector(Ce)) rnorm(n, sd=sqrt(ve*Ce)) else t(chol(ve*Ce)) \%*\% rnorm(n)
Y <- mu + (if (is.vector(G)) G * w else G \%*\% w) + eps

u2 <- system.time(g2 <- genZuchtwert(Y=Y, G=G, Ce = Ce,  Cw=Cw, maxratio=100,
                                   marginal=TRUE))
print(u2)
str(g2)
u <- system.time(g <- genZuchtwert(Y=Y, G=G, Ce = Ce,  Cw=Cw, maxratio=100,
                                 marginal=FALSE))
print(u)
str(g)
str(list(mu=mu, w=w, vw=vw, ve=ve, eps))


hist(g$w - w, 50)
hist(g2$w - w, 50, add=TRUE, border="red")
print(summary(cbind(w, g$w - w, (g$w - w) / w, g2$w - w, (g2$w - w) / w)))

X11()
hist(g$eps - eps, 50)
hist(g2$eps - eps, 50, add=TRUE, border="red")
print(summary(cbind(eps, (g$eps - eps) / eps, (g2$eps - eps) / eps )))

print(rbind(c((g$vw - vw) / vw, (g2$vw-vw)/ vw),
            c((g$ve - ve) / ve, (g2$ve-ve) / ve)))

}

}

\keyword{regression}
