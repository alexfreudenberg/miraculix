\name{codeGeno}
\alias{codeGeno}
\alias{decodeGeno}
\alias{haplo2geno}
\alias{copyGeno}
\alias{zeroNthGeno}
%\alias{decodeNthGeno}
\alias{createSNPmatrix}
\alias{fillSNPmatrix}
\title{
  Creating Compressed Genotype Data
}
\description{
\command{codeGeno} transforms a vector or matrix to a coded matrix of
SNP (genotype) values.
%
\command{decodeGeno} transforms a matrix of coded SNP (genotype) values
into a human readable SNP matrix (or vector).
%
%\command{decodeNthGeno} transforms some selected rows of
%a matrix of coded SNP (genotype) values
%into a human readable SNP matrix (or vector).

%
\command{copyGeno} copies a coded SNP matrix
%
\command{zeroNthGeno} writes zeros into selected rows of coded SNP matrix
%
\command{haplo2geno} transforms the compressed haplotype of one
individuum into its genotype
%
\command{createSNPmatrix} creates an empty compressed \eqn{snps \times
  indiv} matrix.
%
\command{fillSNPmatrix} allows to fill (or replace) colums of a compressed
  \eqn{snps \times indiv} matrix.
}
\usage{
codeGeno(SNPxIndiv)
decodeGeno(SNPxIndiv, do.centering=FALSE, N)
copyGeno(SNPxIndiv)
zeroNthGeno(SNPxIndiv, N)
haplo2geno(CM)
createSNPmatrix(snps, individuals)
fillSNPmatrix(SNPxIndiv, indiv, values)
}
\arguments{
%  \item{G}{\eqn{\{0,1\,2\}}-valued vector or \eqn{snps \times indiv} matrix.  }

  \item{SNPxIndiv}{a compressed SNP (genotype) vector or matrix
    obtained from  \command{codeGeno}}

  \item{do.centering}{
    logical. If \code{TRUE} the value of
    \code{RFoptions()$relationshipmatrix$centering} is considered.
    \describe{
      \item{\code{TRUE}}{
	centering by \command{rowMeans}.
      }
      \item{\code{FALSE}}{
	no centering is performed (although \code{do.centering = TRUE}!)
      }
      \item{\code{is.numeric}}{
	the values given by the user are substracted.
      }
    }
  }
  
  \item{N}{vector of integers, which gives the selected rows. If missing
  all rows are selected.}
  
  \item{CM}{an object returned by \command{\link{codeHaplo}}}

  \item{snps, individuals}{integers, defining the size of the
  compressed \eqn{snps \times indiv} matrix.}

  \item{indiv}{integer vector. It gives
    the column of the \eqn{SNP\times Indiv} matrix
    that has to be filled with \code{values}
  }
  
  \item{values}{coded or uncoded vector or matrix of haplotype or
    genotypes.
  }
}


\details{
  Genotype is coded in 2 bits each locus.

  Technical note: there might be preceding and succeeding bytes
  in the coding to guarantee alignment of the SIMD operations.
}


\value{
  \command{codeGeno} : compressed SNP matrix, i.e. an integer vector with
  attribute \code{"Information"}.

  \command{decodeGeno} : a decompressed i.e., human readable, integer
  vector or \eqn{SNP\times Indiv} matrix. 

  \command{copyGeno} : compressed SNP matrix

  \command{zeroNthGeno} : compressed SNP matrix

\command{haplo2geno} : a genotype vector for one individual

  \command{createSNPmatrix} : a (empty) compressed SNP matrix.

  \command{fillSNPMatrix} : a compressed SNP matrix.
}


\author{\martin}

\seealso{
  \command{\link{vectorGeno}} for multiplying a vector from the left

  \command{\link{genoVector}} for multiplying a vector from the right

}

\examples{ %   library(miraculix)
require(RandomFieldsUtils)

indiv <- sample(1000, 1)
snps <- indiv * 2^sample(7,1)
M <- matrix(nrow = snps, sample(0:2, snps * indiv, replace=TRUE))
storage.mode(M) <- sample(c("integer", "double"), 1)
CM <- codeGeno(M)
str(CM)
Z <- decodeGeno(CM)
Print(M, CM, Z)
stopifnot(all(M == Z))

N <- sample(snps, snps / 4)
Z1 <- decodeGeno(CM, N=N)
stopifnot(all(M[N, ] == Z1))


}

\keyword{misc}
