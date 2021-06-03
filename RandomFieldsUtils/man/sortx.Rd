\name{sortx}
\alias{sortx}
\title{Sorting Vectors}
\description{
  \command{sortx}
  has the same functionality as \command{\link[base]{sort}},
  except that \code{sortx(..., from=from, to=to)}
  is the same as \code{sort[from:to]}

   Sort a vector or factor into ascending or
     descending order. 
}

\usage{
sortx(x, from=1, to=length(x), decreasing=FALSE, na.last = NA) 
}
\arguments{
  \item{x}{an atomic vector}
  \item{from,to}{\code{sort(..., from=from, to=to)} equals
  \code{sort(...)[from:to]}}
  \item{decreasing}{ logical.  Should the sort sort be increasing or
    decreasing? }
  \item{na.last}{for controlling the treatment of \code{NA}s.  If
    \code{TRUE}, missing
    values in the data are put last; if \code{FALSE}, they are put
    first; if \code{NA}, they are removed (see the Notes in
    \command{\link[base]{sort}}) }
}

\value{
 vector of length \code{to}-\code{from}+1.
}
\details{
  The smaller the difference \code{to}-\code{from} is
  compared to the length of \code{x}, the
  faster is \command{sortx} compared to \link[base]{sort}.

  Particularly, \code{sortx(..., from=k, to=k)} 
  is much faster than \code{sort(...)[k]}.
  
  For further details see \link[base]{sort}.
}

\seealso{
  \link{orderx}
}
\examples{
x <- runif(10^6)
k <- 10
system.time(y<-sort(x)[1:k])
system.time(z<-sortx(x, from=1, to=k)) ## much faster
stopifnot(all(y == z)) ## same result
}
\author{Martin Schlather, \email{schlather@math.uni-mannheim.de}
}
\keyword{univar}
\keyword{manip}
