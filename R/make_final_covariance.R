# Authors
# 
#Malena Erbe malena.erbe@agr.uni-goettingen.de 
# Copyright (C) 2013 -- 2014 Malena Erbe
#
# This is proprietary software. All rights belong to the authors.


finalCov<-function(G,A,perc) {
	# following Meuwissen et al. (2011)
	# step 1: calculate Fit_i
	Fit<-diag(G)-1
	# step 2: calculate Fst_G
	FstG<-mean(diag(G))-1
	# step 3: calculate Fst_A
	FstA<-mean(diag(A))-1
	# step 4: calculate Fis_i
	Fis<-(diag(G)-1-FstG)/(1-FstG)
	# step 5: new elements in G
	keepDiag<-diag(G)
	newDiag<-FstA+(1-FstA)*Fis+1
	kin<-(G/2-FstG)/(1-FstG)
	G<-2*(FstA+(1-FstA)*kin)
	diag(G)<-newDiag
	# step 6: combine A and G
	finalCov<-(1-perc/100)*G+(perc/100)*A
	return(finalCov)
}
