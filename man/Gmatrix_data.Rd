\name{Gmatrix_data}
\alias{Gmatrix_data}
\alias{GmatrixData}
\title{
 Gmatrix_data 
}
\description{
  multicore calculation of a genomic relationship matrix following VanRaden 2007,2008 (VanRaden method I) with the genotypes read from an external file
}
\usage{
GmatrixData(zahl=5, filename='file.txt', roundit, nind, nsnps, chr=1,
      asreml=TRUE, freq=NULL, teiler=NULL, teiler.store=NULL, recoded=NULL) 
}
\arguments{
  \item{zahl}{number of threads, default value is 5
  }
  \item{filename}{ name of external input file containing genotypes coded 0,1,2 with SNPs in rows and animals in columns
  }
  \item{roundit}{allele frequencies can be rouded
  }
  \item{nind}{number of individuals
  }
  \item{nsnps}{number of SNPs
  }
  \item{chr}{DOC MISSING
  }
  \item{asreml}{default is TRUE. If TRUE, output format is compatible for ASReml
  }
  \item{freq}{DOC MISSING
  }
  \item{teiler}{DOC MISSING
  }
  \item{teiler.store}{if teiler.store unequals NULL, the value of the denominator is stored
  }
  \item{recoded}{just allowed to be TRUE at the moment 
  }
}

\details{
function uses Fortran code to speed up the calculation of the matrix product ZZ'
 }

\value{
if teiler.store=NULL: G as an R matrix object
if teiler.store unequals NULL: list with Gmatrix in "G" and denominator in "teiler"

}


\author{Malena Erbe, \email{malena.erbe@agr.uni-goettingen.de}
}

\examples{
## to do
}

\keyword{htest}
\keyword{ts}
