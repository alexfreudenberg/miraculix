
ENV <- environment()


## code habpo-type efficiently (for original generation)
codeHaplo <- function(M, N) { ## method checked in C
  H <- if (missing(N)) .Call(C_codeHaplo, M)  else .Call(C_codeHaplo2, M, N)
  attr(H, "method") <- Haplo
  class(H) <- "genomicmatrix"
  H
}
decodeHaplo <- function(CM) {## method checked in C
  method <- attr(CM, "method")
  if (!is(CM, "genomicmatrix") || method !=Haplo)
    stop("matrix not a 'haplo'-typed 'genomicmatrix'")
  .Call(C_decodeHaplo, CM) ## reverse
}


## geno functions
haplo2geno <- function(CM) {## method checked in C
  if (dolocking(NULL))
    stop("if 'dolocking' then implicite re-coding is not possible.")
  if (!is(CM, "genomicmatrix") || attr(CM, "method") != Haplo)
    stop("matrix not a 'haplo'-typed 'genomicmatrix'")
  obj <- .Call(C_haplo2geno, CM)  # CM coded geno matrix
  attr(obj, "method") <- Shuffle
  class(obj) <- "genomicmatrix"
  return(obj)
}

codeGeno <- function(SNPxIndiv) as.genomicmatrix(SNPxIndiv)

decodeGeno <- function(SNPxIndiv, do.centering=FALSE, N) {
  if (!missing(N)) return(decodeNthGeno(SNPxIndiv, N=N, do.centering=FALSE))
  if (is.matrix(SNPxIndiv)) {
    snpXindiv <- SNPxIndiv
    d <- dim(snpXindiv)
    attributes(snpXindiv) <- NULL
    dim(snpXindiv) <- d
  } else {
    stopifnot(is(SNPxIndiv, "genomicmatrix"))
    method <- attr(SNPxIndiv, "method")
    snpXindiv <- (if (!is.null(attr(SNPxIndiv, "filename")))
		    .Call(C_file_get, SNPxIndiv)
		  else if (method  == Shuffle) .Call(C_decodeGeno, SNPxIndiv)
		  else if (method %in% c(TwoBit, ThreeBit))
		    .Call(C_matrix_get, SNPxIndiv)
		  else if (method %in% c(Hamming2, Hamming3))
		    .Call(C_matrixH_get, SNPxIndiv)
		   # stop("transformation to human readable matrix only works for 'Shuffle', 'TwoBit', 'ThreeBit'")
		  else stop("unknown method")
		  )
  }

  if (do.centering) {
    centered <- RFoptions(GETOPTIONS="relationshipmatrix")$centered      
    if (is.na(centered)) return(.Call(C_substract_centered, snpXindiv))
    else if (centered) return(snpXindiv - rowMeans(snpXindiv))
  }
  
  return(snpXindiv)
}

decodeNthGeno <- function(SNPxIndiv, N, do.centering=FALSE) {
  if (any(!is.finite(N) | N < 1))
    stop("all values of 'N' must be positive integers")
  if (is.matrix(SNPxIndiv)) return(SNPxIndiv[N, ])
  N <- as.integer(N - 1)

  stopifnot(is(SNPxIndiv, "genomicmatrix"))
  method <- attr(SNPxIndiv, "method")
  if (method != Shuffle) stop("function only available for method 'Shuffle'")
  
  snpXindiv <-.Call(C_decodeNthGeno, SNPxIndiv, N)

  if (do.centering) {
    centered <- RFoptions(GETOPTIONS="relationshipmatrix")$centered      
    if (is.na(centered)) return(.Call(C_substract_centered, snpXindiv))
    else if (centered) return(snpXindiv - rowMeans(snpXindiv))
  }
  
  return(snpXindiv)
}

copyGeno <-function(SNPxIndiv) {
  if (is.matrix(SNPxIndiv)) return(SNPxIndiv + 0)
  stopifnot(is(SNPxIndiv, "genomicmatrix"))
  method <- attr(SNPxIndiv, "method")
  if (method != Shuffle) stop("function only available for method 'Shuffle'")
  obj <- .Call(C_copyGeno, SNPxIndiv)
  attr(obj, "method") <- Shuffle
  class(obj) <- "genomicmatrix"
  return(obj)
}

zeroNthGeno <- function(SNPxIndiv, N) {
  if (any(!is.finite(N) | N < 1))
    stop("all values of 'N' must be positive integers")
  if (is.matrix(SNPxIndiv)) {
    SNPxIndiv[N, ] <- 0
    return(SNPxIndiv)
  }
  N <- as.integer(N - 1)

  stopifnot(is(SNPxIndiv, "genomicmatrix"))
  method <- attr(SNPxIndiv, "method")
  if (method != Shuffle) stop("function only available for method 'Shuffle'")

  snpXindiv <- copyGeno(SNPxIndiv)
#  str(snpXindiv)
  .Call(C_zeroNthGeno, snpXindiv , N)
  snpXindiv 
}

createSNPmatrix <- function(snps, individuals) {## method checked in C
  ## note instead of snps, also a human readable matrix can be provided
  if (dolocking(NULL))
    stop("if 'dolocking' then implicite coding is not possible.")
  if (is.matrix(snps)) {
    stopifnot(missing(individuals))
    C <- .Call(C_codeGeno, snps)
  } else C <-.Call(C_createSNPmatrix, as.integer(snps), as.integer(individuals))
  attr(C, "method") <- Shuffle 
  class(C) <- "genomicmatrix"
  C
}


fillSNPmatrix <- function(SNPxIndiv, indiv, values) {## method checked in C
  ## SNPxIndiv[, indiv] <- values
  ## note: values can be anything: haplo, geno, coded or uncoded
  if (!is(values, "genomicmatrix")) {
    if (is.matrix(values) && nrow(values) == 2)
      values <- codeHaplo(values)
    else if (is.vector(values) || (is.matrix(values) && ncol(values) == 1))
      values <- .Call(C_codeGeno, values)
    else stop("unknown kind of 'values'")
  } 

   Infos <- attr(values, "information")
  p <- Infos[[INFO_INFO+1]]
 
  if (length(p) > 0 && p[1] == HAPLO){
    values <- haplo2geno(values)
  }
  else if (p[1] != GENO) stop("strange coded values")
  .Call(C_fillSNPmatrix, SNPxIndiv, as.integer(indiv), values)
}


relationshipMatrix <- function(SNPxIndiv, ...) {
  if (length(list(...)) > 0) {
    RFoptOld <- RFoptions(## LOCAL=1,
                          SAVEOPTIONS="relationshipmatrix", ...)
    on.exit(RFoptions(LIST=RFoptOld))
  }
  RFopt <- RFoptions(GETOPTIONS = c("basic", "relationshipmatrix"))

  centered <- RFopt$relationship$centered
  opt <- RFopt$basic
  

  if (!is(SNPxIndiv, "genomicmatrix")) SNPxIndiv <- as.genomicmatrix(SNPxIndiv)
  Infos <- attr(SNPxIndiv, "information")
  method <- attr(SNPxIndiv, "method")

  if (opt$verbose) cat("method =", RELSHIP_METH_NAME[method + 1],
		       "using", opt$cores, "core(s).\n")

#  str(cSNPxIndiv)
#  Print(Infos)

  if (method == Shuffle) 
    A <- .Call(C_relationshipMatrix, SNPxIndiv)

  else if (method==TwoBit || method == ThreeBit) {
##    print(SNPxIndiv)
    if (RFopt$relationship$per_snp) {
      compressed <- Infos[[INFO_INFO + 1]][MEM + 1]
      ## C coding of indices:
      snp.slices <- as.integer(round(compressed / opt$cores * 0:(opt$cores-1)))
      snp.slices2 <- c(snp.slices[-1], compressed)
      stopifnot(is.integer(snp.slices2)) 
      
      ## MmPsqSNPs needs hash table
      parfun <- function(i)
	.Call(C_MmPsqSNPs, SNPxIndiv, snp.slices[i], snp.slices2[i])
      z <- parallel::mclapply(1:opt$cores,  mc.cores=opt$cores, parfun)
      if (!is.list(z) || any(sapply(z, class) == "try-error")) stop(z)
      A <- .Call(C_getRelmatrixSNP, z, SNPxIndiv)
    } else {
      A <- .Call(C_MmPsqIndividualsParallel, SNPxIndiv)
    }

  } else if (method %in% c(Hamming2, Hamming3))
    A <- .Call(C_matrix_mult, SNPxIndiv)    

  else if (method == NoSNPcoding) {
    p <- rowMeans(SNPxIndiv)
    if (!is.na(centered) && centered) SNPxIndiv <- SNPxIndiv - p
    if (RFopt$relationship$normalized)
      SNPxIndiv <- SNPxIndiv / sqrt(sum(p * (1 - p / 2)))
    A <- base::crossprod(SNPxIndiv)
  }

  else stop("method '", method, "' not programmed yet.")
  
  if (is.na(centered)) {
##    print("isna centered")
    p  <- .Call(C_get_centered)
    pA <- vectorGeno(p, SNPxIndiv=SNPxIndiv, do.centering=FALSE)
    A <- A - pA - rep(1, length(pA)) %*% t(pA) + sum(p^2)
  }
  attr(A, "centered") <- centered
  return(A)
}


vectorGeno <- function(V, SNPxIndiv, do.centering=FALSE) {
  if (do.centering) stop("centered not programmed yet")
  ## RFopt <- RFoptions(GETOPTIONS="relationshipmatrix")
  if (is.matrix(SNPxIndiv)) return(V %*% SNPxIndiv)
  stopifnot(is(SNPxIndiv, "genomicmatrix"))
  method <- attr(SNPxIndiv, "method")
  if (method == Shuffle) .Call(C_vectorGeno, V, SNPxIndiv)
  else if (any(method == c(TwoBit, ThreeBit)) &&
	   !is.null(attr(SNPxIndiv, "filename"))) {    
###   if (is.matrix(SNPxIndiv)) {
###    snpXindiv <- g %*% SNPxIndiv   
###    if (is.na(centered))
###      snpXindiv <- snpX indiv - sum(.C(C_get_centered) * g)
###    else if (centered) snpXindiv <- snpXindiv - sum(rowMeans(SNPxIndiv) * g)
###    return(snpXindiv)
####  }
    return(.Call(C_dot_file, SNPxIndiv, as.double(V)))
  }

  as.vector(V %*% as.matrix(SNPxIndiv))
}


genoVector <- function(SNPxIndiv, V, do.centering=FALSE) {
  if (do.centering) stop("centered not programmed yet")
  if (is.matrix(SNPxIndiv)) return(SNPxIndiv * V)
  stopifnot(is(SNPxIndiv, "genomicmatrix"))
  method <- attr(SNPxIndiv, "method")
  if (method == Shuffle) return(.Call(C_genoVector, SNPxIndiv, V))
  else if (any(method == c(TwoBit, ThreeBit)) &&
	   !is.null(attr(SNPxIndiv, "filename"))) {    
   return(.Call(C_file_dot, SNPxIndiv, as.double(V)))
  } 

  as.matrix(SNPxIndiv) %*% V
}


allele_freq <- function(SNPxIndiv) {
  # name <- as.character(as.list(match.call(envir=parent.frame(2L)))[[1]])
  if (is.character(SNPxIndiv)) SNPxIndiv <- as.genomicmatrix(SNPxIndiv)
  if (!is(SNPxIndiv, "genomicmatrix")) {
    if (!is.matrix(SNPxIndiv))
      stop("not a matrix and not of class 'genomicmatrix'")
    if (nrow(SNPxIndiv) < ncol(SNPxIndiv))
      message("are you sure you have give an snp x individual matrix?")
    return(rowMeans(SNPxIndiv) * 0.5)
  }
  method <- attr(SNPxIndiv, "method")
  if (method == Shuffle) return(.Call(C_allele_freq, SNPxIndiv))

  rowMeans(as.matrix(SNPxIndiv)) * 0.5
}

solveRelMat <- function(A, tau, vec, betahat, destroy_A=FALSE) {
  ## if not destroyed, more memory is consumed
  ## implements ( t(Z) Z + tau \1_{ind x ind} )^{-1} %*% vec
  stopifnot(length(tau) == 1 && length(betahat) ==1)
  .Call(C_solveRelMat, A, as.double(tau), as.double(vec), as.double(betahat),
	destroy_A)
}

                
SNPeffect <- function(SNPxIndiv, vec, centered=TRUE, tau=0) {
  snpXind <- SNPxIndiv
  geno <- is(snpXind, "genomicmatrix")
  if (!geno) snpXind <- as.genomicmatrix(snpXind)
  rM <- relationshipMatrix(snpXind, centered=centered, normalized=FALSE)
  if (tau != 0) rM <- rM + diag(nrow(rM)) * tau ## ehemals eps
  rM <- solvex(rM, vec)

  method <- attr(SNPxIndiv, "method")
  return(if (geno) SNPxIndiv %*% rM else genoVector(SNPxIndiv, rM))
}

