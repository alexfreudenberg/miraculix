
# Authors
#
# Martin Schlather, schlather@math.uni-mannheim.de
# Copyright (C) 2013 -- 2015  Martin Schlather
#               2014 Florian Skene
#

as.matrix.genomicmatrix <- function(x, ...) decodeGeno(x, ...)

print.genomicmatrix <- function(x, ...) {
  method <- attr(x, "method")
  cat("genomic matrix (method = ", RELSHIP_METH_NAME[method + 1], "): ",
      sep="")
  Infos <- attr(x, "information")
  info <- Infos[[INFO_INFO + 1]]
  if (is.null(fn <- attr(x, "filename"))) {
    cat(sep="", " int [1:", info[SNPS + 1], ", 1:", info[INDIVIDUALS + 1],
	  "]\n")
  } else {
    cat(sep="", " int [1:", info[SNPS + 1], ", 1:", info[INDIVIDUALS + 1],
	  "] <", fn, ">\n")
  }
}

str.genomicmatrix <- function(object, ...) {
  print.genomicmatrix(object, ...)
  cat("Attributes are a ")
  str(attributes(object)) ##
}

as.genomicmatrix <- function(obj, obj.type,
			     coding, header, individualspercolumn,
			     doubledindividuals, leadingcolumns, ...) {

  if (length(list(...)) > 0) {
    RFoptOld <- RFoptions(##LOCAL=1,
                          SAVEOPTIONS="relationshipmatrix", ...)
    on.exit(RFoptions(LIST=RFoptOld))
  }
  RFopt <- RFoptions(GETOPTIONS=c("basic", "relationshipmatrix"))
  opt <- RFopt$relationship
  basic <- RFopt$basic
  method <- opt$relship_method
  if (method == AutoCoding)
    if (.Call(C_sse_check)) method <- Shuffle
    else {
      RFoptions(per_snp = FALSE)
      method <- ThreeBit
    }
  if (is(obj, "genomicmatrix")) {
    warning("'obj' is already of class 'genomicmatrix'. Its original definition is kept.");
    return(obj)
  }

  if (is.numeric(obj)) {
    if (nrow(obj) <= ncol(obj))
      warning("The number of SNPs ist less than or equal to the number of individuals")
    if (!is.integer(obj)) {
##      Print(obj)
      storage.mode(obj) <- "integer"
      if (basic$helpinfo)
	message("Better use 'integer' as storage mode for 'SNPxIndiv") 
    }
    if (any(obj != 0  && obj != 1 && obj != 2)) stop("not a genomic matrix")
    if (method == Shuffle) {
      if (dolocking(NULL))
	stop("if 'dolocking' then implicite coding is not possible.")
      obj <- .Call(C_codeGeno, obj)
    } else if (method == TwoBit || method == ThreeBit) {
      obj <- .Call(C_matrix_coding, obj)
    } else if (any(method == c(Hamming2, Hamming3))) {
     ## factor <- if (method == Hamming2) 96 else 128
      ##if ((mis <- nrow(obj %% factor)) != 0) {
	  ##	obj <- rbind(obj, matrix(0L, factor - mis, ncol(obj)))
	  ##	stopifnot(is.integer(obj))
	  ## }
      ## str(obj)
      obj <- .Call(C_matrix_codingH, obj)
    } else if (method == NoSNPcoding)  {
      attr(obj, "information") = list(
	info = c(
	  what = NA,
	  snps = nrow(obj),
	  individuals=ncol(obj),
	  addr0=0,
	  addr1=0,
	  addr2=0,
	  addr3=0,
	  N = NA,
	  isSNPxInd = TRUE,
	  blocks = NA,
	  bits = NA,
	  bitsperblock = NA,
	  snpspercompressed = NA,
	  header=NA,
	  doubledindividuals=NA,
	  leadingcolumns=NA
	),
	code = NULL,
	p = NULL,
	pmi = NULL,
	coding = NULL
      )
    } else stop("unknown method")
    
    attr(obj, "method") <- method
    class(obj) <- "genomicmatrix"
    return(obj)
  }

  if (!is.character(obj)) stop("'SNPxIndiv' must be a file name or a matrix")
   
  individuals <- snps <- NA
  if (missing(obj.type)) {
    split <- strsplit(obj, "\\.")[[1]]
    file.ext <-  c("txt", "bgl", "phased", "tped", "ped", "bed")
    s <- pmatch(split[length(split)], file.ext)
    if (!is.na(s)) { 
      obj.type <- c(txt="beagle", bgl="beagle", phased="plink",
                     tped="plink2", ped="plink2",
                     bed="plinkbinary")[file.ext[s]]
    } else {
      message("Note: obj.type could not be detected.")
    }
    }
  
  Coding <-c(beagle="AB? ", plink="AB? ", plink2="12? ", plinkbinary="12345")
  Header <-    c(beagle=1,    plink=0,    plink2=0,     plinkbinary=-3)
  AniPerCol <- c(beagle=TRUE, plink=TRUE, plink2=FALSE, plinkbinary=TRUE)
  DoubledAni <-c(beagle=TRUE, plink=TRUE, plink2=FALSE, plinkbinary=FALSE)
  LeadingCol <-c(beagle=2,    plink=4,    plink2=6,     plinkbinary=0)
  
  if (!missing(obj.type)) {
    if (missing(coding)) coding <- Coding[obj.type]
    if (missing(header)) header <- Header[obj.type]
    if (missing(individualspercolumn)) individualspercolumn <- AniPerCol[obj.type]
    if (missing(doubledindividuals)) doubledindividuals <- DoubledAni[obj.type]
    if (missing(leadingcolumns)) leadingcolumns <- LeadingCol[obj.type]
    
    if (obj.type == "plinkbinary") {
      s <- paste(split[-length(split)], collapse=".")
      
      file.lines <- function(obj) {
        sys <- strsplit(system(paste("wc", obj), intern=TRUE), " ")[[1]]
        sys <- sys[nchar(sys) > 0]
          as.integer(sys[1])
      }
      
      ##  #lines of the files give the snps or individuals
      individuals <- file.lines(paste(s, ".fam", sep=""))
      snps <- file.lines(paste(s, ".bim", sep=""))
    }
  }
  if (nchar(coding) != 4 && nchar(coding) != 5)
    stop("haplotype must be given as 4 characters, e.g. 'AB? '")
  if (nchar(coding) == 4) {
    ## never change this coding
    coding <- 
      paste(substr(coding, 1, 2), substr(coding, 2, 4),
            sep="")
  }
  
  if (missing(coding) || missing(header) || missing(individualspercolumn) ||
      missing(doubledindividuals) || missing(leadingcolumns)) {
    if (missing(coding) && missing(header) && missing(individualspercolumn) &&
        missing(doubledindividuals) && missing(leadingcolumns))
      stop("unknown obj.type. Please give the obj.type")
    stopifnot(!missing(coding))
    stopifnot(!missing(header))
    stopifnot(!missing(individualspercolumn))
    stopifnot(!missing(doubledindividuals))
    stopifnot(!missing(leadingcolumns))
  }

  obj <- filename <- as.character(obj) 
  attr(obj, "information") = list( ## this list must match exactly INFO_*
    info = c(
      what = NA,
      snps = as.integer(snps),
      individuals=as.integer(individuals),
      addr0=0,
      addr1=0,
      addr2=0,
      addr3=0,
      N = NA,
      isSNPxInd = individualspercolumn,
      blocks = NA,
      bits = NA, ## 10
      bitsperblock = NA,
      snpspercompressed = NA,
      header=as.integer(header),
      doubledindividuals=as.integer(doubledindividuals),
      leadingcolumns=as.integer(leadingcolumns)
    ),
    code = NULL,
    p = NULL,
    pmt = NULL,
    coding = as.character(coding)
  )

  obj <- .Call(c(C_file_Shuffle, C_file_coding, C_file_coding, C_file_Hamming,
		 C_file_Hamming, C_failure, C_file_get, C_failure)[method + 1],
	       obj)
  attr(obj, "method") <- method
  attr(obj, "filename") <- filename
  class(obj) <- "genomicmatrix"
  return(obj)
}

getMT <- decodeGeno
MTdot <- genoVector
dotMT <- vectorGeno
relM <- relationshipMatrix
snpEffect <- SNPeffect

crossprodx <- function(SNPxIndiv) relationshipMatrix(SNPxIndiv, centered =FALSE)
