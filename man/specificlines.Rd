\name{SpecificLines}
\alias{SpecificLines}
\title{
 SpecificLines
}
\description{
reads only specific lines (defined by user) of a file 
}
\usage{
SpecificLines(infile,ntotal,whichlines,idcols,datacols)
}
\arguments{
  \item{infile}{name of input file
  }
  \item{ntotal}{number of total lines of the file
  }
  \item{whichlines}{ numeric vector stating the lines that should be read 
  }
  \item{idcols}{ number of leading columns in the input file containing additional info 
  }
  \item{datacols}{ number of columns containing data in the input file
  }
}

\details{
uses a Fortran code to speed up the reading  
}

\value{ 
returns a list with values for leading info columns in "id" and for data columns in "datamat" 
}


\author{Malena Erbe, \email{malena.erbe@agr.uni-goettingen.de}
}

\examples{

# to be done

}

\keyword{IO}
