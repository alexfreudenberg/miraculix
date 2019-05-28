

# Authors
# 
#Malena Erbe malena.erbe@agr.uni-goettingen.de 
# Copyright (C) 2013 -- 2015 Malena Erbe
#
# This is proprietary software. All rights belong to the authors.



Gmatrix<-function(zahl=5,genos.gen,roundit=3,asreml=TRUE,freq=NULL,teiler=NULL) {
  if (is.null(freq)) {
    cat('berechne Frequenzen','\n')
    freq<-round(apply(genos.gen,2,mean)/2,roundit)
  }
  cat('erstelle Z','\n')
  Z<-sweep(genos.gen,2,2.0 * freq,'-')
  #library(multicore)  ## martin:sollte ohne gehen
  nsnps<-dim(Z)[2]
  nind<-dim(Z)[1]

                                                                                 ##zahl<-5
  cat('berechne Einzelteile von G','\n')
  teile<-mclapply(1:zahl,function(x){
    nn <- as.integer(floor(nsnps/zahl))
    y <- double((nind+1)*nind/2)
     if (x<zahl) {
      cat(x,'\n')
      xx<-t(Z[,((x-1)*nn+1):(x*nn)])
      .Fortran(C_gmatrix, xx, y, nn, nind)
    } else if (x==zahl) {
      cat(x,'\n')
      xx<-t(Z[,((x-1)*nn+1):nsnps])
      .Fortran(C_gmatrix, xx, y, as.integer(nsnps-nn*(zahl-1)), nind)
    }
  }
                  ,mc.cores=zahl)
  cat('verbinde Einzelteile','\n')
  for (j in 1:zahl) {
    if (j==1) loesung<-teile[[j]][[2]] ## geaendert von Malena Mitte Jan 2014
    if (j!=1) loesung<-loesung+teile[[j]][[2]]
  }
  G<-matrix(0,nind,nind)
  G[lower.tri(G,diag=TRUE)]<-loesung
  G<-G+t(G)-diag(diag(G))
  if (is.null(teiler)) {
    G<-G/sum(2*freq*(1-freq))
  } else {
    G<-G/teiler
  }
  if (asreml) {
    cat('bereite G fuer ASReml vor','\n')
    Glong<-G[upper.tri(G,diag=TRUE)]
    Glong<-cbind(rep(1:nind,1:nind),sequence(1:nind),Glong)
    return(Glong)
  } else {
    return(G)
  }
}
