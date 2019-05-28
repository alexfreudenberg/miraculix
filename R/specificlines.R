# Authors
# 
#Malena Erbe malena.erbe@agr.uni-goettingen.de 
# Copyright (C) 2013 -- 2014 Malena Erbe
#
# This is proprietary software. All rights belong to the authors.


SpecificLines<-function(infile,ntotal,whichlines,idcols,datacols) {
  stop("function must be rewritten")
  ntotal <- as.integer(ntotal)
  whichlines <- as.integer(whichlines)
  idcols <- as.integer(idcols)
  datacols <- as.integer(datacols)
	
  write.table(infile,'paramfile_specificlines',col.names=FALSE,
	      row.names=FALSE,quote=FALSE,sep='\t')

  nlines<-length(whichlines)
  datamat<-matrix(0L,nrow=nlines,ncol=datacols)
  
  dataresults<-.Fortran(C_specific_lines, ntotal, nlines, whichlines,
			datacols, idcols, datamat)		

	idresults<-read.table('idfile_specificlines',header=FALSE,colClasses='character')
	
	results<-list(datamat=dataresults[[6]],id=idresults)
	
	system('rm idfile_specificlines paramfile_specificlines')
	
	return(results)
}
