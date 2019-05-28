\name{vectorGeno}
\alias{vectorGeno}
\alias{genoVector}
\title{
  Multiplication of a vector to a compressed SNP matrix
}
\description{
\command{vectorGeno} multiplies a vector from the left onto a compressed
SNP matrix.

\command{genoVector} does it from the right.
}
\usage{
vectorGeno(V, SNPxIndiv, do.centering=FALSE)
genoVector(SNPxIndiv, V, do.centering=FALSE)
}
\arguments{
  \item{SNPxIndiv}{a compressed SNP (genotype) vector or matrix
 obtained from  \command{codeGeno}. Uncoded SNP matrix is also possible.}

 \item{V}{numerical vector}

 \item{do.centering}{not programmed yet.
   %see \command{\link{decodeGeno} }
 }
}


\details{
  Let $G$ be the (original) \eqn{SNP\times Indiv} matrix.
  \command{vectorGeno} and \command{genoVector} return \eqn{VG} and
  \eqn{GV}, respectively.
}

\value{
  vector of length \code{nrow(CG_SNPxIndiv)} and
  \code{ncol(CG_SNPxIndiv)}
  for \command{vectorGeno} and \command{genoVector}, respectively.
}


\author{\martin}

\examples{ %   library(miraculix)
require(RandomFieldsUtils)

indiv <- 1 + sample(500, 1)
snps <- indiv * 2^sample(7, 1)
snps <- indiv * 100
M <- matrix(ncol=indiv, sample(0:2, indiv * snps, replace=TRUE))
print(system.time(CM <- codeGeno(M)))


## V %*% G
Vl <- runif(snps)
print(system.time(VM1 <- vectorGeno(Vl, CM))) ## slower than standard
print(system.time(VM <- as.vector(Vl \%*\% M)))
stopifnot(all.equal(VM, VM1))

## G %*% V
Vr <- runif(indiv)
print(system.time(MV1 <- genoVector(CM, Vr))) ## faster than standard
print(system.time(MV <- as.vector(M \%*\% Vr)))
stopifnot(all.equal(MV, MV1))
}

\keyword{misc}
