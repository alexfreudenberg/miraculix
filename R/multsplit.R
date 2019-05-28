# Authors
# 
#Malena Erbe malena.erbe@agr.uni-goettingen.de 
# Copyright (C) 2013 -- 2014 Malena Erbe
#
# This is proprietary software. All rights belong to the authors.

multSplit<-function(eins,zwei,zahl) {

	n<-dim(eins)[1]
	nn<-floor(n/zahl)
	teile<-mclapply(1:zahl,function(x){
		if (x<zahl) {
#			cat(x,'\n')
			xx<-eins[((x-1)*nn+1):(x*nn),]
		} else if (x==zahl) {
#			cat(x,'\n')
			xx<-eins[((x-1)*nn+1):n,]    
		}
		xx%*%zwei
	}
	)
	
	for (j in 1:zahl) {
		if (j==1) loesung<-teile[[j]]
		if (j!=1) loesung<-rbind(loesung,teile[[j]])
	}
	
	return(loesung)
	
}
