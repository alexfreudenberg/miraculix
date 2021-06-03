\name{RFoptions}
\alias{RFoptions}
\title{Setting control arguments}
\description{
 \command{\link{RFoptions}} sets and returns control arguments for
 diverse packages (\pkg{miraculix}, \pkg{RandomFields}).

 \command{\link{RFoptions}} should not be used within parallelizing R
 commands such as \command{mclapply} in package \pkg{parallel}.
% \command{\link[parallel]{mclapply}}
}
%\usage{
% RFoptions(..., no.readonly = TRUE)
%}
%\arguments{
% \item{...}{arguments in \code{tag = value} form, or a list of tagged values.}
% \item{no.readonly}{If \command{\link{RFoptions}} is called without
%   argument then all arguments are returned in a list. If
%   \code{no.readonly=TRUE} then only rewritable arguments are returned.
%   Currenlty all arguments are rewritable. So the list is empty.
% }
%}
\details{
%  Note that this function can only be used after
%  \code{library(RandomFieldsUtils)} has been made available explicitely.

  The specific parameters for \pkg{miraculix} are the following.
  \RFU in \pkg{RandomFieldsUtils} for further options.
  
  \describe{
    \item{efficient}{ logical.
      If \code{TRUE} then always the most time efficient code is used.

      Default: \code{TRUE}.
      It is strongly recommended to retain this value. 
    }
      
    \item{centered}{logical or numerical.
      If \code{TRUE} the \eqn{P} matrix is substracted before the
      crossproduct of the the SNP matrix is calculated,
      see \code{\link{relationshipMatrix}} for the \eqn{P} matrix.

      If numeric, then the length of this vector must equal the number
      of SNPs per individual. Then this vector is substracted for each
      individual. 
      Furthermore, \code{normalized} is FALSE.
      As the size of \code{centered} can be large, this vector is never
      returned by \code{RFoption()}; instead \code{NA} is returned.

      Note that \code{centered} also sets the value of \code{normalized}.
      
      Default : \code{TRUE}
    }
   

    \item{\code{cores}}{
      Number of cores for multicore algorithms.
      
      Default : 1
    }

   \item{digits}{OBSOLETE. scalar.
      If \code{digits} is negative no rounding is performed.
      Else the matrix \eqn{P} when calculating the relationsship matrix
      \eqn{(M-P)^T(M-P)} is rounded to the given number of absolute (not
      significant) digits.
      
      Default : \code{3.0}.
    }

    \item{mma_method}{integer.

      Default : \code{4}. It is strongly recommended to retain this value.  
    }
 
    \item{normalized}{logical.
      If \code{TRUE} the relationship matrix is normalized by
      \eqn{\sigma^2}, see \code{\link{relationshipMatrix}}.
      

      Its value is set to the value of \code{centered} whenever
      the value of \code{centered} is changed. So \code{normalized}
      must be set always after \code{centered},
      e.g. \code{RFoptions(centered=TRUE, normalized=FALSE)}, but not
      \code{RFoptions(normalized=FALSE, centered=TRUE)}.
      
      Default : \code{TRUE}
    }

  
    \item{snpcoding}{integer. This parameter influences memory
      requirements and the speed of the algorithms.

      Default : \code{TwoBit}.
      All man pages assume that the default value has been retained.
      So, it should not changed.
    }

    
    \item{returnsigma}{
      logical. Whether \eqn{\sigma^2} shall be also returned when
      the relationship matrix is calculated.
    }
  }  
}






\value{
 \code{NULL} if any argument is given, and the full list of
 arguments, otherwise.
}

%\references{}
\author{Martin Schlather, \email{schlather@math.uni-mannheim.de}
 \url{http://ms.math.uni-mannheim.de/de/publications/software}}

\seealso{
  \link{benchmarking} for the effect of changing
  the values of \code{efficient}, \code{mma_method} and \code{snpcoding}.
  
  \RFU
  
}


\examples{% library(adoption)
RFoptions()$genetics
}

\keyword{models}

