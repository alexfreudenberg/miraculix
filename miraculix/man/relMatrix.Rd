\name{relationshipMatrix}
\alias{relationshipMatrix}
\alias{solveRelMat}
\alias{allele_freq}
\alias{vwm}
\alias{SNPeffect}
\alias{crossprodx}


\title{
  Fast calculation of the Genomic Relationship Matrix
  and its derivatives
}
\description{
\command{relationshipMatrix} calculates the relationship matrix
\eqn{A=(M-P)^T (M-P) /\sigma^2} from the SNP matrix \eqn{M}
where \eqn{P = p (1,\ldots,1)} with  \eqn{p = M \%*\% (1,\dots,1)^T / n}.
 Furthermore, \eqn{sigma^2} equals
  \eqn{\sigma^2 = p^T  (1 - p/2)\in[0,\infty)}.


\command{crossprodx} calculates the cross-product of \code{SNPxIndiv},
i.e. it is identical to call \command{relationshipMatrix} with
optional argument,
\code{centered=FALSE}, cf. \command{\link{RFoptions}}


\command{allele_freq} calculates \eqn{p/2}.

\command{SNPeffect} calculates \eqn{M (A + \tau I)^{-1} v}

\command{solveRelMat} calculates
\deqn{(A + \tau I)^{-1} v} and \deqn{A(A + \tau I)^{-1} v + \beta}
where \eqn{A} is the relationship matrix.

}

\usage{
relationshipMatrix(SNPxIndiv, ...)
crossprodx(SNPxIndiv) 

solveRelMat(A, tau, vec, betahat=NULL, destroy_A=FALSE)
SNPeffect(SNPxIndiv, vec, centered=TRUE, tau=0)
allele_freq(SNPxIndiv)
}

\arguments{
  \item{SNPxIndiv}{\eqn{\{0,1\,2\}}{0,1,2}-valued (snps \eqn{\times} indiv)
    matrix or the result of \command{\link{genomicmatrix}}.
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
  \item{betahat}{scalar or \code{NULL}. See also section value.}
   \item{destroy_A}{logical. If \code{TRUE} the values of the
    matrix \code{A} will be overwritten
    during the calculations (leading to a faster execution with
    less memory needs).}
}


\details{
  Let \eqn{p = M \%*\% (1,\dots,1)^T / n} where \eqn{n} is the number
  of individuals.
  Then, the matrix \eqn{P} equals \eqn{P = p (1,\ldots,1)}.

  The constant \eqn{sigma^2} equals
  \eqn{\sigma^2 = p^T  (1 - p/2)}.
  
  \code{solveRelMat} has a speed and memory advantage in comparison to the
  direct implementation of the above formulae.
}



\value{
  \command{relationsshipMatrix} returns a (Indiv \eqn{\times} Indiv)
  numerical matrix.

  The return value of \command{solveRelMat} depends on \code{betahat}.
  If the latter is \code{NULL}, only the vector \eqn{(A + \tau I)^{-1} v}
  is returned. Else,  a list of 2 elements is returned. First element equals
  the vector \deqn{(A + \tau I)^{-1} v,} the
  second element equals \deqn{A(A + \tau I)^{-1} v + \beta.}
}


\me

\seealso{
  \code{\link{RFoptions}} for options,

  \link{benchmarking} for implementation details.
   
  
}

\examples{ %   library(miraculix)
set.seed(0)

indiv <- 1 + sample(100:500, 1)
snps <- indiv * 2^sample(1:if (interactive()) 7 else 5, 1)
M <- matrix(ncol=indiv, sample(0:2, indiv * snps, replace=TRUE))
alleleFreq <- 0.5 * rowMeans(M)
print(system.time(G <- genomicmatrix(M)))
print(G)  

## crossprodx vs crossprod: about 10x faster
Print(system.time(C <- crossprodx(M)))   
print(system.time(C2 <- crossprod(M)))
Print(C, C2, range(C - C2))
stopifnot(all(C == C2))

## relationshipMatrix vs crossprod: about 10x faster
Print(system.time(R <- relationshipMatrix(M)))
print(system.time(R <- relationshipMatrix(G)))
print(system.time({
  sigma2 <- 2 * sum(alleleFreq * (1 - alleleFreq))
  R2 <- crossprod(M - 2 * alleleFreq) / sigma2
}))
stopifnot(all.equal(as.double(R), as.double(R2)))

}

\keyword{misc}
