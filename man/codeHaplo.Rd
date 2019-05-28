\name{codeHaplo}
\alias{codeHaplo}
\alias{decodeHaplo}
\title{
  Compressing Haplotype Data
}
\description{
\command{codeHaplo} codes
}
\usage{
codeHaplo(M, N)
decodeHaplo(CM)
}
\arguments{
  \item{M, N}{\eqn{\{0,1\}}-valued.
    Either \code{M} is a \eqn{2\times SNP}{2 x SNP} matrix coding the
    two alleles of an individual. Or \code{M} and \code{N} are the
    two vectors for the two alleles.
  }
  \item{CM}{a vector obtained from coding a matrix by \command{codeHaplo}
  }
}


\details{
  The bitwise coding is (locus 1, allele 1), (locus 1, allele 2), (locus 2,
  allele 1) etc.

  Technical note: there might be preceding and succeeding bytes
  in the coding to guarantee alignment of the SIMD operations.
}


\value{
  \command{codeHaplo} : a vector with length equal to the number of
  rows of \code{M}.

  \command{decodeHaplo} : a decompressed, i.e. human readable integer vector of 4 components.
}


\author{\martin}

\examples{ %   library(miraculix)
require(RandomFieldsUtils)

n <- sample(1000, 1)
M <- matrix(nrow = 2, sample(0:1, 2 * n, replace=TRUE))
storage.mode(M) <- sample(c("integer", "double", "logical"), 1)
CM <- codeHaplo(M)
Z <- decodeHaplo(CM)
Print(M, CM, Z)
stopifnot(all(M == Z))

}

\keyword{misc}
