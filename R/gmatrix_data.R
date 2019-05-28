

# Authors
# 
# Malena Erbe malena.erbe@agr.uni-goettingen.de 
# Copyright (C) 2013 -- 2014 Malena Erbe
#
# This is proprietary software. All rights belong to the authors.


GmatrixData<-function(zahl=5, filename='file.txt', roundit, nind, nsnps, chr=1, asreml=TRUE, freq=NULL, teiler=NULL, teiler.store=NULL, recoded=NULL) {
 
  nn<-as.integer(floor(nsnps/zahl))
  nind <- as.integer(nind)
  snps <- as.integer(snps)
  roundit <- as.integer(roundit)
	
  cat('berechne Einzelteile von G','\n')
  cat(nn,'\n')
  teile<-mclapply(1:zahl, function(x){	
    xx<-paste(chr,'000',x,sep='')
    write.table(matrix(c(filename,x),nrow=1),paste('param',xx,sep=''),col.names=FALSE,row.names=FALSE,quote=FALSE,sep='\t')
    y <- double((nind+1)*nind/2)
    runnr<-as.integer(xx)
    nenner<-double(1)
    cat(x,'\n')
    paramfile <- paste("param", runnr, sep="")
    if (is.null(recoded)) {
      ##cat('oder hier','\n')
      .Fortran(C_gmatrix_data, nsnps, nind, nn, runnr, roundit, y, nenner,
	       paramfile)
    } else {	
      ##cat('hallo hier','\n')
      .Fortran(C_gmatrix_data_recoded, nsnps, nind, nn, runnr, roundit,
	       y, nenner, paramfile)
    }
  }
, mc.cores=zahl)
  
  cat('verbinde Einzelteile','\n')
  for (j in 1:zahl) {
    if (j==1) {
      loesung<-teile[[j]][[6]]
      teiler<-teile[[j]][[7]]
    } else if (j!=1) {
      loesung<-loesung+teile[[j]][[6]]
      teiler<-teiler+teile[[j]][[7]]
    }
  }
  ## return(teile)
  G<-matrix(0,nind,nind)
  G[lower.tri(G,diag=TRUE)]<-loesung
    G<-G+t(G)-diag(diag(G))
  if (is.null(teiler.store)) {
    G<-G/teiler
	} else {
	  G<-list(G=G,teiler=teiler)
	}

  if (asreml) {
    if (!is.null(teiler.store)) {
      stop('asreml=TRUE and teiler.store=TRUE together do not make sense!')
    }
    cat('bereite G fuer ASReml vor','\n')
    Glong<-G[upper.tri(G,diag=TRUE)]
    Glong<-cbind(rep(1:nind,1:nind),sequence(1:nind),Glong)
    return(Glong)
  } else {
    return(G)
  }
}
