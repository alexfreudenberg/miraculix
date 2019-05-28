\name{FileMatrix}
\alias{FileMatrix}
\title{
 Fast calculation of \eqn{(M-P) (M-P)^T /\sigma^2}
}
\description{
  This function calculates the relationship matrix
}\usage{
FileMatrix(file, file.type,  coding, header,
    animalspercolumn, doubledanimals, leadingcolumns,
    centered = TRUE)
}
\arguments{
  \item{file}{Character for a file, or a matrix.
    
    The file format depends
    on the extension as well as on
    \code{file.type}, \code{coding}, \code{header},
    \code{animalspercolumn}, and \code{leadingcolumns}

    If \code{file.type} is missing, the \code{file.type} is
    \itemize{
      \item \code{"beagle"} if the file extension is in \code{"txt",
	"bgl", "phased"}. In this case we have\cr
      \code{coding = "AB? "}, \code{header=1}, \code{animalspercolumn=TRUE}, 
      \code{doubledanimals=TRUE}, \code{leadingcolumns=2}
      
      \item\code{"plink"} if the file extension equals \code{"tped"}.
      In this case we have\cr
      \code{coding = "AB? "}, \code{header=0}, \code{animalspercolumn=TRUE}, 
      \code{doubledanimals=TRUE}, \code{leadingcolumns=4}
      
      \item\code{"plink2"} if the file extension equals \code{"ped"}.
      In this case we have\cr
      \code{coding = "12? "}, \code{header=0}, \code{animalspercolumn=FALSE}, 
      \code{doubledanimals=FALSE}, \code{leadingcolumns=6}
    }
  }

   \item{file.type}{character: \code{"beagle", "plink"}.
    if not \code{missing} and the respective arguments are missing the
    following values are set
    \itemize{
      \item\code{"beagle"}
      \itemize{
	\item \code{coding="AB? "}
	\item \code{header=1}
	\item \code{animalspercolumn=TRUE}
	\item \code{leadingcolumns=2}
      }
      \item\code{"plink"}
      \itemize{
	\item \code{coding="AB? "}
	\item \code{header=0}
	\item \code{animalspercolumn=TRUE}
	\item \code{leadingcolumns=4}
      }
    }
  }
  \item{coding}{
    string with 4 or 5 characters with the following meaning
    \itemize{
      \item first character: haplotype or genotype coded with \code{0}
      \item second character: haplotype or genotype coded with \code{1}
      \item middle character: if given than a genotype matrix is
      indicated and the middle character defines \code{2}
      \item last but one character: missing value;
      if there are missing values the algorithm currently fails.      
      \item last character: space character
    }
    if any other character is found in the genuine columns of the
    matrix, the function fails.

    See also \code{animalpercolumn}.
   }
  \item{header}{
    number of header lines
  }
  \item{animalspercolumn}{logical.
    If \code{TRUE} a SNPs x INDIVIDUAL matrix is expected in case
    of genotypes and SNPs x (2 * INDIVIDUALS) matrix in case of
    haplotypes.
    
    If \code{FALSE} an INDIVIDUUM x SNPs matrix is expected in case
    of genotypes and (2 * INDIVIDUALS) x SNPs matrix in case of haplotypes.
  }
  \item{doubledanimals}{logical.
    if \code{TRUE} and haplotypes are given, then two rows (or columns)
    of animals define the haplotype. Otherwise two columns (or rows)
    of SNPS define the haplotype.
  }
  \item{leadingcolumns}{
    number of first columns to be ignored
  }
  \item{centered}{
    if \code{TRUE} or \code{is.numeric}
    then \eqn{M-P} instead of Matrix \eqn{M} is
    considered.

    If \code{is.numeric} then \code{centered} must have length
    equal to the number of snps.
    }
}

%\details{}
 
\value{
  It returns a \eqn{m\times m} matrix where \eqn{m} equals the number
  of individuals involved.
}

\details{
  The matrix \eqn{P} equals
  \eqn{P = (1,\ldots,1)^T  colMeans(M)}
  where \eqn{M^T} equals \code{MT}.

  The constant \eqn{sigma^2} equals
  \eqn{\sigma^2 = 2 P_{1, \cdot} (1- P_{1,\cdot})^T}.
  
}

\author{Martin Schlather, \email{schlather@math.uni-mannheim.de}
}

\examples{

## see relM

}

\keyword{htest}
\keyword{ts}
