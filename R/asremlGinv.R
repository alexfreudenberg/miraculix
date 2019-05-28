## Authors
## 
##Malena Erbe malena.erbe@agr.uni-goettingen.de 
## Copyright (C) 2015 Malena Erbe
##
## This is proprietary software. All rights belong to the authors.
##
## solve changed to solvex 6.6.18 by MS


asremlGinv<-function(G,variable='ID',rown=NULL) {

	Gi<-solvex(G)
	nn<-dim(Gi)[1]
	Ginv<-data.frame(row=rep(1:nn,1:nn),column=sequence(1:nn),x=Gi[upper.tri(Gi,diag=TRUE)])
	colnames(Ginv)[3]<-variable
	if (is.null(rown)) {
		attr(Ginv,'rowNames')<-1:nn
	} else {
		attr(Ginv,'rowNames')<-rown
	}
	return(Ginv)	
}

