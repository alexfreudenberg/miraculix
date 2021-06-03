## This script benchmarks miraculix against R-base in different environments. CPU and GPU versions are tested. We simulate random SNP matrices.
# In the second part we benchmark the empirical decomposition of phenotypic variances as laid out by Schreck (2019)
library(miraculix)
library(MoBPS)

# Set global options
RFoptions(cores = 8, efficient=FALSE)
indiv_sizes <- as.integer(seq(5e3,50e3, length.out = 10))
snps <- 20e3
# Matrix for time measurements
times <- matrix(0,ncol=length(indiv_sizes), nrow=3)
colnames(times) <- as.character(indiv_sizes)
rownames(times) <- c("R","Shuffle","MMAGPU")

# iterate over number of individuals
for (indiv in indiv_sizes){

  # Simulate SNP matrix
  M <- matrix(ncol=indiv, rep(sample(0:2, 1e3 * snps, replace=TRUE), indiv/1e3))
  storage.mode(M) <- "integer"
  cat("\n\n")
  
  # Start time benchmarks
  print(times[1,as.character(indiv)] <- system.time(crossprod(M))[3])
  RFoptions(snpcoding = Shuffle256)
  
  print(times[2,as.character(indiv)] <- system.time(crossprodx(M))[3])
  RFoptions(snpcoding = MMAGPU)
  {crossprodx(M[1:1e3,1:1e3]);NULL}
  print(times[3,as.character(indiv)] <- system.time(crossprodx(M))[3])
  
  gc()
  
}
# Remove buffers and save results
rm(M)
save.image("benchmarks_GPU.RData")

# A shorter time benchmark is provided in the following. We use the median of 5 iterations to get stable results. 
num_cols <- c(2e2,1e3,5e3,10e3,20e3)
# Prepare time matrix
speed_gpu <- speed_avx <- matrix(ncol = length(num_cols),nrow=5, dimnames = list(NULL,num_cols))
RFoptions(snpcoding=MMAGPU,cores=12, la_mode=LA_GPU )
for(n in num_cols){
  snps <- n * 4
  # Simulation
  M <- matrix(ncol=n, sample(0:2, n * snps, replace=TRUE))
  G <- genomicmatrix(M)
  # Warm-up run
  relationshipMatrix(G)
  # Benchmark
  for(i in 1:5)speed_gpu[i,match(n,colnames(speed))] <- system.time(relationshipMatrix(G))[3]
}

RFoptions(snpcoding=Shuffle256, la_mode=LA_INTERN)
for(n in num_cols){
  snps <- n * 4
  # Simulation
  M <- matrix(ncol=n, sample(0:2, n * snps, replace=TRUE))
  G <- genomicmatrix(M)
  # Warm-up run
  relationshipMatrix(G)
  # Benchmark
  for(i in 1:5)speed_avx[i,match(n,colnames(speed))] <- system.time(relationshipMatrix(G))[3]
}


## Plot function
df_speed <- data.frame(index = num_cols, gpu = apply(speed_gpu,2,median), avx = apply(speed_avx,2,median))
ggplot(tidyr::pivot_longer(df_speed,cols=c("gpu","avx")), aes(x = index, y = log10(value), colour = name))+geom_line() + scale_y_continuous(n.breaks = 10)


# Benchmark decomposition of phenotypic variance 
indiv_sizes <- as.integer(seq(1e3,20e3))

times <- matrix(0,ncol=length(indiv_sizes), nrow=2)
colnames(times) <- as.character(indiv_sizes)
rownames(times) <- c("sommer","miraculix")
nsnp <- 10000
nenv <- 500
source("~/master_thesis/decomp_schreck.R")
# Iterate
for(n in indiv_sizes){
  # Simualate population
  population <- creating.diploid(nsnp = nsnp, nindi = n, n.additive = c(100,100))
  population <- breeding.diploid(population, heritability = c(0.5,0.5),
                                 phenotyping.database = cbind(1,1),
                                 n.observation = c(1,0))
  
  population <- breeding.diploid(population, heritability = c(0.5,0.5),
                                 phenotyping.database = cbind(1,2),
                                 n.observation = c(1,1))
  
  population <- breeding.diploid(population, bve=TRUE,
                                 bve.gen = 1)
  
  
  Z <- get.geno(population, gen=1) # Bestimmung der Genotypen
  pheno <- get.pheno(population, gen=1) # Bestimmung der Phaenotypen
  genomische_wert <- get.bv(population, gen=1) # Bestimmung tatsaechlicher genomischer Zuchtwerte (diese sind in der Praxis unbekannt)
  G <- miraculix::relationshipMatrix( miraculix::genomicmatrix(Z), centered = TRUE, normalized = TRUE) # Berechnung Verwandtschaftsmatrix
  X <- matrix(rnorm(nenv * n),nrow=n)
  beta <- matrix(rnorm(nenv),ncol=1)
  y <- matrix(pheno[1,],ncol=1)+ X%*%beta
  # Start time benchmarks
  times[1,as.character(n)] <-   system.time(print(Var_decomp_lmm(y,X,t(Z))$total))[3]  
  times[2,as.character(n)] <- system.time(print(cu_decomp(y,X,t(Z))$total))[3]
  
}
