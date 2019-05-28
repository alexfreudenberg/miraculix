\name{Gmatrix}
\alias{Gmatrix}
\title{
 Gmatrix 
}
\description{
  multicore calculation of a genomic relationship matrix following VanRaden 2007,2008 (VanRaden method I)
}
\usage{
Gmatrix(zahl=5,genos.gen,roundit=3,asreml=TRUE,freq=NULL,teiler=NULL) 
}
\arguments{
  \item{zahl}{number of threads, default value is 5
  }
  \item{genos.gen}{marker matrix of dimension individuals x SNPs, genotypes coded 0/1/2, no missing values 
  }
  \item{roundit}{allele frequencies can be rouded, default value is 3 
  }
  \item{asreml}{gives output in format compatible for ASReml
  }
  \item{freq}{if vector of allele frequencies is given here, these frequencies are used for correction and scaling
  }
  \item{teiler}{if a numeric value is given here, this value will be used for standardisation 
  }
}

\details{to do 
 }

\value{to do 

}


\author{Malena Erbe, \email{malena.erbe@agr.uni-goettingen.de}
}

\examples{
## to do
}

\keyword{htest}
\keyword{ts}
