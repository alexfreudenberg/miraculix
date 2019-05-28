\name{finalCov}
\alias{finalCov}
\title{
 finalCov
}
\description{
  combination of G and A matrices based on paper of Meuwissen et al. (2011) (J Anim Breed Genet. 2011 Dec;128(6):429-39. doi: 10.1111/j.1439-0388.2011.00966.x)
}
\usage{
finalCov(G,A,perc)
}
\arguments{
  \item{G}{genomic relationship matrix in R matrix format
  }
  \item{A}{pedigree relationship matrix in R matrix format 
  }
  \item{perc}{numeric value between 0 and 100 being the weight in percentage that should be given to A 
  }
}

\details{
	The following steps are performed:
	1. calculate F_{it}(i)
	2. calculate F_{st}(G)
	3. calculate F_{st}(A)
	4. calculate F_{is}(i)
	5. build updated G matrix (G_{new})
	6. combine A and new G matrix as (1-perc/100)*G_{new}+(perc/100)*A
 }

\value{ returns the final covariance matrix in R matrix format as a combination of A and updated G

}


\author{Malena Erbe, \email{malena.erbe@agr.uni-goettingen.de}
}

\examples{
if (require(synbreedData) && require(pedigree)) {
data(cattle)

# genotypes
genos<-cattle$geno
genos[genos=='AA']<-0
genos[genos=='AB']<-1
genos[genos=='BB']<-2

storage.mode(genos)<-'integer'

# allele frequencies
freq<-round(apply(genos,2,mean,na.rm=TRUE)/2,3)

# Z matrix
Z<-sweep(genos,2,2*freq,'-')
Z[is.na(Z)]<-0

# build G matrix
G<-tcrossprod(Z)/sum(2*freq*(1-freq))

# pedigree
ped<-cattle$pedigree

# build A matrix
makeA(ped[,1:3],is.element(ped$ID,rownames(G))) # take care: A.txt is written!
Amat<-matrix(scan('A.txt'),ncol=3,byrow=TRUE)
A<-matrix(0,dim(G)[1],dim(G)[1])
A[cbind(rep(1:dim(G)[1],1:dim(G)[1]),sequence(1:dim(G)[1]))]<-Amat[,3]
A<-A+t(A)-diag(diag(A))

# order A according to G
orderA<-match(rownames(G),ped$ID[is.element(ped$ID,rownames(G))])
A<-A[orderA,orderA]

# combine G and A with 10% A and 90% G
GA<-finalCov(G,A,10)
}
}

\keyword{htest}
\keyword{ts}
