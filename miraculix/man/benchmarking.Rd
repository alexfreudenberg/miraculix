\name{Benchmarking}
\alias{benchmarking}
\alias{Shuffle}
\alias{Shuffle256}
\alias{TwoBit}
\alias{Packed}
\alias{Packed256}
\alias{Multiply}
\alias{Multiply256}
\alias{ThreeBit}
\alias{Hamming2}
\alias{Hamming3}
\alias{AutoCoding}
\alias{NoSNPcoding}
\alias{NoSNPcodingR}
\alias{FirstGenuineMethod}
\alias{FirstMoBPSmethod}
\alias{LastGenuineMethod}
\alias{LastMoBPSmethod}
\alias{SNPCODING_NAMES}
\title{Benchmarking for Calculation the Genomic Relationship Matrix}

\description{
  When the standard options
  (\code{\link{RFoptions}()$genetics$efficient=TRUE} and 
  \code{\link{RFoptions}()$genetics$snpcoding=TwoBit}) are retained,
  \pkg{miraculix} uses a compressed
  2-bit representation of the genomic values and detects the most
  efficient available code according to the hardware and the compiler
  flags. (See the installation hints at \link{miraculix}.

  This man page is for users who are interested in coding details of
  miraculix. In particular, we give benchmarks for (essentially)
  all the implemented variants.
}

\details{
  \describe{
    \item{efficient}{ logical.
      If \code{TRUE} then always the most time efficient code
      is used for the given \code{snpcoding}, independently what is specified
      by the user. For instance,
      \itemize{
	\item{\code{TwoBit} (no SIMD needed)}
	\item{\code{Packed} (SSE2 needed)}
	\item{\code{Shuffle} (SSSE3 needed)}
	\item{\code{Shuffle256} (AVX2 needed)}
	\item{\code{MMAGPU} (GPU (CUDA needed)}
      }
      are based on the same SNP coding. So these 5 options are
      indistinguishable if \code{efficient=TRUE}.

      The value \code{FALSE} is mainly for mixed architectures
      and benchmarking.
    }
 
      
    \item{snpcoding}{integer. 
      Possible values are
      \describe{
	\item{\code{MMAGPU}}{GPU (CUDA)}
	\item{\code{Shuffle}}{two bit mini hash table based on SSSE3}
	\item{\code{Shuffle256}}{two bit mini hash table based on AVX2}
	\item{\code{Packed}}{4-bit integer arithmetic based on SSE2}
	\item{\code{Packed256}}{4-bit integer arithmetic based on AVX2}
	\item{\code{Multiply}}{16-bit integer arithmetic based on SSE2}
	\item{\code{Multiply256}}{16-bit integer arithmetic based on AVX2}
	\item{\code{Twobit}}{two bit hash table}
	\item{\code{Threebit}}{three bit hash table}
	\item{\code{Hamming2}}{4-bit representation; used in PLINK}
 	\item{\code{Hamming3}}{4-bit representation; used in PLINK}
	\item{\code{AutoCoding}}{method is chosen by the programme itself}
	\item{\code{NoSNPcoding}}{no coding, i.e. 32 bit integer. It
	  used AVX2, if available}
	\item{\code{NoSNPcodingR}}{No coding: 32 bit integer, R
	  code. Only for testing purposes.}
      }

      In \code{for} loops that run through all available methods
      the constants \code{FirstGenuineMethod} and
      \code{LastGenuineMethod}
      might be useful, see below.
      (In case of the package \pkg{MoPBS} or if interest is in the 2 bit methods
      only, use the constants \code{FirstMoBPSmethod} and
      \code{LastMoBPSmethod}.)

      In case the names of the method is needed, use 
      \code{SNPCODING_NAMES[snp_coding + 1]}.      
     }
   }
 }

 \note{The installation for both SIMD and GPU must be performed
  manually with the CRAN version of miraculix; see the directory
  inst/GPU or GPU for the 
  necessary files to be exchanged (do not forget to run autoconf)
  }

 
\section{Benchmarks}{
  Computing times for the relationship matrix in comparison to 'crossprod'
  in standard implementation on a NVIDIA Graphic Card % alex !! TO DO
  and an Intel(R) Core(TM) i7-8550U CPU @
  1.80GHz, R version 3.6.0 (Linux) with \code{indiv = 1000} and
  \code{snps = 5e5} are:
  MMAGPU        : >100 faster (CUDA; 16x compressed; \code{mma_method=4})\cr 
  MMAGPU        : >100 faster (CUDA; 8x compressed on GPU; \code{mma_method=2})\cr 
  MMAGPU        : >100 faster (CUDA; x compressed; \code{mma_method=1})\cr 
  Shuffle256    : 48 x faster (AVX2; 16x compressed)\cr
  Packed256     : 36 x faster (AVX2; 16x compressed)\cr
  Shuffle       : 35 x faster (SSSE3; 16x compressed)\cr
  Multiply256   : 29 x faster (AVX2; 16x compressed)\cr
  Packed        : 28 x faster (SSE2; 16x compressed)\cr
  Hamming2      : 24 x faster (SSE2; 4x compressed)\cr
  Hamming3      : 21 x faster (SSSE3; 4x compressed)\cr
  Multiply      : 17 x faster (SSE2; 16x compressed)\cr
  ThreeBit      : 17 x faster (uint64_t; 10x compressed)\cr
  TwoBit        : 15 x faster (uint64_t; 16x compressed)\cr
  NoSNPcoding   :  4 x faster (int, AVX2; not compressed)\cr
  NoSNPcodingR  : calls \link[base]{crossprod} 
 
 In parantheses,
 first the instruction set or s the main data type is given,
 then the data compression with respect to 32 bit integer.
 
The following code was used: % library(miraculix)
\preformatted{
RFoptions(cores = 1, efficient=FALSE)
indiv <- 1000
snps <- 5e5 ## may cause memory allocation problems in R; better use 5e4 !!
methods <- c(NoSNPcodingR, FirstGenuineMethod:LastGenuineMethod)
M <- matrix(ncol=indiv, sample(0:2, indiv * snps, replace=TRUE))
for (storageMode in c("integer", "double")){
    storage.mode(M) <- storageMode
  cat("\n\n")
  print(S <- system.time(C <- crossprod(M)))
  p <- rowMeans(M)
  P <- p \%*\% t(rep(1, indiv))
  sigma2 <- sum(p * (1- p/2))
  A <- crossprod(M-P) / sigma2
  print(S <- system.time(C <- crossprod(M)))  
  for (method in methods) {
    RFoptions(snpcoding = method)
    cat("\nstorage=", storageMode, "  method=", SNPCODING_NAMES[method + 1],
    "\n")
    S0 <- system.time(G <- genomicmatrix(M))
    print(S1 <- system.time(C1 <- crossprodx(M)))
    print(S2 <- system.time(C2 <- crossprodx(G)))
    stopifnot(all(C == C1))
    stopifnot(all(C == C2))
    R1 <- S / S1 
    R2 <- S / S2 
    print(0.5 * (R1 + R2))
    print(S3 <- system.time(C3 <- relationshipMatrix(M)))
    print(S4 <- system.time(C4 <- relationshipMatrix(G)))
    R3 <- S / S3
    R4 <- S / S4 
    print(0.5 * (R3 + R4))
    stopifnot(all.equal(as.double(A), as.double(C3)))
    stopifnot(all.equal(as.double(A), as.double(C4)))
    gc()
  }
}
}
}

\examples{ %   library(miraculix)
require(RandomFieldsUtils)
set.seed(0)
snpcodes <- c(TwoBit, ThreeBit)
if (has.instruction.set("SSE2")) snpcodes <- c(snpcodes, Hamming2)
if (has.instruction.set("SSSE3")) snpcodes <- c(snpcodes, Hamming3, Shuffle)
if (has.instruction.set("AVX2")) snpcodes <- c(snpcodes, Shuffle256)
if (has.instruction.set("CUDA")) snpcodes <- c(snpcodes, MMAGPU)
   
Print(snpcodes)

indiv <- 1 + sample(100:500, 1)
snps <- indiv * 2^sample(1:if (interactive()) 7 else 5, 1)
RFoptions(snpcoding=sample(snpcodes, 1))
M <- matrix(ncol=indiv, sample(0:2, indiv * snps, replace=TRUE))
print(system.time(G <- genomicmatrix(M)))
print(G)  

## crossprodx vs crossprod: about 10x faster
Print(system.time(C <- crossprodx(M)))
print(system.time(C2 <- crossprod(M)))
stopifnot(all(C == C2))

## allele_freq vs rowMeans: about equally fast
Print(system.time(af <- allele_freq(M)))
print(system.time(alleleFreq <- 0.5 * rowMeans(M)))
stopifnot(all.equal(as.double(alleleFreq), as.double(af)))

## relationshipMatrix vs crossprod: about 10x faster
Print(system.time(R <- relationshipMatrix(M)))
print(system.time(R <- relationshipMatrix(G)))
print(system.time({
  sigma2 <- 2 * sum(alleleFreq * (1 - alleleFreq))
  R2 <- crossprod(M - 2 * alleleFreq) / sigma2
}))
stopifnot(all.equal(as.double(R), as.double(R2)))

### solveRelMat vs. solve: about equally fast
tau <- 0.0001
vec <- runif(indiv)
beta <- runif(1)
Print(system.time(S <- solveRelMat(R, tau=tau, vec=vec, betahat=beta)))
print(system.time({r <- solve(R + diag(indiv) * tau, vec);
                   y <- as.vector(R \%*\% r + beta)}))
stopifnot(all.equal(S$rest, r))
stopifnot(all.equal(S$yhat, y))

}

\keyword{algebra}
