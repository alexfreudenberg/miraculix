## This script illustrates further examples for benchmarking miraculix and RandomFieldsUtils
# Should be understood as a collection of code snippets

suppressMessages({
library(miraculix)
require(RandomFieldsUtils)
set.seed(5L)
})
#dyn.load("sample_snpmatrix.so")
RFoptions(cores=8,max_chol=30e3,useGPU=T) 
indiv <- 10e3
snps <- indiv*4
cat("Generating SNP matrix\n")
#print(system.time({M <- .C("samplesnps", as.integer(indiv * snps), as.integer(seed), M = integer(indiv * snps))$M;M <- matrix(ncol=indiv, M)}))
M<-matrix(sample(0:2,indiv *snps, replace=T),ncol=indiv)
cat(paste0("Size of M: ", format(object.size(M),units = "Gb"), "\n"))

RFoptions(centered=F, normalized=F)
RFoptions(snpcoding=NoSNPcodingR)
##Warmup
# cat("Warmup\n")
# suppressMessages({relationshipMatrix(matrix(ncol=50,nrow=20, sample(0:2,1000,replace=T)));NULL})

# cat("Encoding:\n")
# print(system.time(G <- genomicmatrix(M)))
# cat("Cublas:\n")
# print(system.time(R2 <- relationshipMatrix(G)))

RFoptions(snpcoding=Shuffle256)
cat("Encoding:\n")
print(system.time(G <- genomicmatrix(M)))
cat("Shuffle:\n")
print(system.time(R <- relationshipMatrix(G)))
# print(all.equal(as.double(R), as.double(R2),tolerance=1e-10,scale=1))



RFoptions(snpcoding=MMAGPU)
reps <- 5
times<-double(reps)

for(i in 1:(reps+1)){
    print(reps[i]<-system.time({
    cat("Encoding:\n")
    print(system.time(G <- genomicmatrix(M)))
    cat("MMAGPU:\n")
    print(system.time(R2 <- relationshipMatrix(G)))
    stopifnot(all.equal(as.double(R), as.double(R2),tolerance=1e-10,scale=1))
}))
}

# for(i in 1:5){
# tau <- 0.0001
# vec <- runif(indiv)
# beta <- runif(1)
# print(system.time(S <- solveRelMat(R2, tau=tau, vec=vec, betahat=beta)))
# }
# print(system.time({r <- solve(R2 + diag(indiv) * tau, vec);y <- as.vector(R2 %*% r + beta)}))
# cat(paste('Deviation of ',sum(abs(S$rest -r))+ sum(abs(S$yhat - y)),'\n'))
# stopifnot(all.equal(S$rest, r))
# stopifnot(all.equal(S$yhat, y))
