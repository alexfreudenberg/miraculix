
## code "all" information of the animals effectively
codeOrigins <- function(M) .Call(C_codeOrigins, M)
decodeOrigins <- function(CM, row)  .Call(C_decodeOrigins, CM, row)



computeSNPS <- function(population, gen, sex, nr,
			from_p = 1, to_p = Inf,
			output_compressed=FALSE,			
			select = NULL,
			what = c("geno", "haplo")
			) {
  ##  Print(gen, sex, nr)
##  str(population$breeding[[gen[1]]][[sex[1]]][[nr[1]]])
##  name <- as.character(as.list(match.call(envir=parent.frame(2L)))[[1]])
##  DEBUG(name, population, gen, sex, nr);
  
  if (is.matrix(gen)) {
    gen <- gen[, 1]
    sex <- gen[,2]
    nr <- gen[, 3]
  }
  geno <- what[1] != "haplo"
 
  ans <-  .Call(C_computeSNPS, population,
		as.integer(gen),
		as.integer(sex),
		as.integer(nr), as.integer(from_p),
		as.integer(if (is.finite(to_p)) to_p else NA), select,
		geno)
  attr(ans, "method") <- if (geno) Shuffle else Haplo
  class(ans) <- "genomicmatrix"
  if (!output_compressed) {
    ##print(ans);    Print(output_compressed)
    
    ans <- if (geno) decodeGeno(ans) else decodeHaplo(ans)
  }
  
  ##  if (PB(name)) RETURN(name) else SAVE(name, ans)
  return(ans)
}


compute <- function(population, gen, sex, nr, tau, vec, betahat,
		    select = NULL, matrix.return=FALSE) {
  ## implements ( t(Z) Z + tau \1_{ind x ind} )^{-1} %*% vec

    
  if (is.matrix(gen)) {
    gen <- gen[, 1]
    sex <- gen[,2]
    nr <- gen[, 3]
  }
  stopifnot(length(tau) == 1 && length(betahat) == 1)
  
   .Call(C_compute, population,
	as.integer(gen),
	as.integer(sex),
	as.integer(nr),
	as.double(tau),
	as.double(vec),
	as.double(betahat),
	as.integer(select),
	as.logical(matrix.return)
	)
}

init_parallelcompute <- function() .Call(C_init_parallelcompute)
