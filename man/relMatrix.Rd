\name{relationshipMatrix}
\alias{relationshipMatrix}
\alias{solveRelMat}
\alias{allele_freq}
\alias{relM}
\alias{vwm}
\alias{getMT}
\alias{MTdot}
\alias{dotMT}
\alias{SNPeffect}
\alias{crossprodx}
\alias{as.matrix.genomicmatrix}
\alias{print.genomicmatrix}
\alias{str.genomicmatrix}


\title{
  Fast calculation of the Genomic Relationship Matrix
}
\description{
\command{relationshipMatrix} calculates the relationship matrix
\eqn{A=(M-P)^T (M-P) /\sigma^2} from the SNP matrix \eqn{M}
where \eqn{P = p (1,\ldots,1)} with  \eqn{p = M \%*\% (1,\dots,1)^T / n}.
 Furthermore, \eqn{sigma^2} equals
  \eqn{\sigma^2 = p^T  (1 - p/2)\in[0,\infty)}.

\command{allele_freq} calculates \eqn{p/2}.

\command{SNPeffect} calculates \eqn{M (A + \tau I)^{-1} v}

\command{solveRelMat} calculates
\deqn{(A + \tau I)^{-1} v} and \deqn{A(A + \tau I)^{-1} v + \beta}
where \eqn{A} is the relationship matrix.

\command{crossprodx} calculates the cross-product of \code{SNPxIndiv}
}

\usage{
relationshipMatrix(SNPxIndiv, ...)
crossprodx(SNPxIndiv) 

solveRelMat(A, tau, vec, betahat, destroy_A=FALSE)
SNPeffect(SNPxIndiv, vec, centered=TRUE, tau=0)
allele_freq(SNPxIndiv)
}

\arguments{
  \item{SNPxIndiv}{\eqn{\{0,1\,2\}}-valued  \eqn{snps \times indiv}
    matrix or the result of \command{\link{codeGeno}}.
  }
  \item{...}{see \command{\link{RFoptions}} -- better use
  \command{\link{RFoptions}}. 
    The main two options are:

    \code{centered}: see below
  
    \code{normalized}:logical. if \code{FALSE} then the division by
    \eqn{sigma^2} is not performed
  }
  \item{centered}{if \code{FALSE} then \eqn{P} is not
    substracted.}
  \item{A}{a symmetric, positive definite matrix, which is a
    relationship matrix}  
  \item{tau}{non-negative scalar}
  \item{vec}{the vector \eqn{v}}
  \item{betahat}{a scalar}
  \item{destroy_A}{logical. If \code{TRUE} the matrix \code{A} will be
    destroyed during the calculations (leading to a faster execution with
    less memory needs).}
}


\details{
  Let \eqn{p = M \%*\% (1,\dots,1)^T / n} where \eqn{n} is the number
  of individuals.
  Then, the matrix \eqn{P} equals \eqn{P = p (1,\ldots,1)}.

  The constant \eqn{sigma^2} equals
  \eqn{\sigma^2 = p^T  (1 - p/2)\in[0,\infty)}.

  \code{solveRelMat} has a speed and memory advantage in comparison to the
  direct implementation of the above formulae.
}



\value{
  \command{relationsshipMatrix} returns \eqn{Indiv \times Indiv}
  numerical matrix.

  \command{solveRelMat} returns a list of 2 elements (depending on
  the value of \code{matrix}. First element equals \deqn{(A + \tau I)^{-1}
    v}. Second element equals \deqn{A(A + \tau I)^{-1} v + \beta}. 
}


\author{\martin}

\examples{ %   library(miraculix)
require(RandomFieldsUtils)

indiv <- 1 + sample(100:500, 1)
snps <- indiv * 2^sample(1:7, 1)
indiv <- 500
snps <- indiv * 2^7
M <- matrix(ncol=indiv, sample(0:2, indiv * snps, replace=TRUE))
storage.mode(M) <- sample(c("integer", "double"), 1)
print(system.time(R <- relationshipMatrix(M)))
print(system.time(af <- allele_freq(M)))
print(system.time({
  alleleFreq <- 0.5 * rowMeans(M) ## p / 2
  sigma2 = 2 * sum(alleleFreq * (1 - alleleFreq))
  R2 <- crossprod(M - 2 * alleleFreq) / sigma2
}))
stopifnot(all.equal(as.double(R), as.double(R2)))
stopifnot(all.equal(as.double(alleleFreq), as.double(af)))


tau <- 0.0001
vec <- runif(indiv)
beta <- runif(1)
print(system.time(S <- solveRelMat(R, tau=tau, vec=vec, betahat=beta)))
print(system.time({r <- solve(R + diag(indiv) * tau, vec);
                   y <- as.vector(R \%*\% r + beta)}))
stopifnot(all.equal(S$rest, r))
stopifnot(all.equal(S$yhat, y))
}

\keyword{misc}
