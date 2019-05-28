\name{multSplit}
\alias{multSplit}
\title{
 multSplit
}
\description{
multicore matrix multiplication by splitting the first matrix into n parts
}
\usage{
multSplit(eins,zwei,zahl)
}
\arguments{
  \item{eins}{first matrix in R matrix format
  }
  \item{zwei}{second matrix in R matrix format 
  }
  \item{zahl}{number of cores to be used (= number of parts the first matrix is splitted into 
  }
}

\details{
uses mcapply from package parallel to perform a multicore matrix multiplication 
}

\value{ 
returns the solution of the matrix multiplication in R matrix format 
}


\author{Malena Erbe, \email{malena.erbe@agr.uni-goettingen.de}
}

\examples{

# two matrices
A<-matrix(1:200,ncol=40)
B<-matrix(320:1,nrow=40)

# normal matrix multiplication
AB <- A\%*\%B

# multicore version on five cores
AB2 <- multSplit(A, B, 5)

# check if results are equal
# Print(AB, AB2)
stopifnot(all(AB == AB2)) 

}

\keyword{htest}
\keyword{ts}
