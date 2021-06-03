## This script benchmarks miraculix against rrblup and base R. CPU and GPU versions are tested. We simulate populations with 50k SNPs and three generations.
library(rgpu)
library(MoBPS)
library(miraculix)

# Number of individuals
n <- c(1000,2000,5000,10000,15000,20000, 30e3,40e3) 
# Prevent timeouts from R
setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE)
# Global RFutils options
RFoptions(cores=8,helpinfo=FALSE)
# Matrix for time measurements
comp_time <- matrix(0, nrow=length(n), ncol=4)
rownames(comp_time) <- n
colnames(comp_time) <- c("rrBLUP","miraclux_CPU","base_R","miraculix_GPU")

#Iterate over number of individuals
for(i in 1:length(n)){

   #Simulate population
  population <- creating.diploid(nsnp = 50e3, nindi = n[i], n.additive = c(100,100))

  population <- breeding.diploid(population, heritability = c(0.5,0.5),
                                 phenotyping.database = cbind(1,1),
                                 n.observation = c(1,0))
  population <- breeding.diploid(population, heritability = c(0.5,0.5),
                                 phenotyping.database = cbind(1,2),
                                 n.observation = c(1,1))

  # Extract relevant information
  Z <- get.geno(population, gen=1) # Bestimmung der Genotypen
  pheno <- get.pheno(population, gen=1) # Bestimmung der Phaenotypen
  genomische_wert <- get.bv(population, gen=1) # Bestimmung tatsaechlicher genomischer Zuchtwerte (diese sind in der Praxis unbekannt)
 
  y <- pheno[1,] # Betrachte Trait 1
  y_real <- genomische_wert[1,]
  ratio <- 1 # Verhaeltnis zwischen residualer/zufaelliger Effekte und genetischer Varianz

  # Start time measurements
  cat("rrBLUP\n")
  RFoptions(snpcoding=NoSNPcodingR)

  comp_time[i,1] <-system.time({
    G <- miraculix::relationshipMatrix( miraculix::genomicmatrix(Z), centered = TRUE, normalized = TRUE) # Berechnung Verwandtschaftsmatrix
    rrBLUP::mixed.solve(y-mean(y), K = G, method="REML", bounds = c(1e-9,1e9))
  })[3]
  
  cat("miraculix useGPU=F\n")
  RFoptions(snpcoding=Shuffle256, la_mode=LA_INTERN)

  comp_time[i,2] <-system.time({
    G <- miraculix::relationshipMatrix( miraculix::genomicmatrix(Z), centered = TRUE, normalized = TRUE) # Berechnung Verwandtschaftsmatrix
    y_hat <- miraculix::solveRelMat(G, ratio, y - mean(y), mean(y))
    })[3] 

  cat("base R\n")
  RFoptions(snpcoding=NoSNPcodingR)

  comp_time[i,3] <-system.time({
    G <- miraculix::relationshipMatrix( miraculix::genomicmatrix(Z), centered = TRUE, normalized = TRUE) # Berechnung Verwandtschaftsmatrix
   ( G %*% (chol2inv(chol(add.diag(G, ratio))) %*% (y-mean(y)))) + mean(y)
  })[3]
  
  cat("miraculix useGPU=T\n")
  RFoptions(la_mode=LA_GPU,snpcoding=MMAGPU)

  comp_time[i,4] <- system.time({
    G <- miraculix::relationshipMatrix( miraculix::genomicmatrix(Z), centered = TRUE, normalized = TRUE) # Berechnung Verwandtschaftsmatrix
    rm(Z); rm(population); gc()
    y_hat2 <- ifelse(n[i]<30e3, miraculix::solveRelMat(G, ratio, y - mean(y), mean(y)), rgpu::cu_matmul(G, rgpu::solve_large(G + ratio *diag(rep(1,n[i])), matrix(y-mean(y),ncol=1) )) + mean(y)  )
  })[3] 

  # Sanity check results and print intermediary results
  print(paste0("Deviation of ",sum(abs(y_hat$yhat-y_hat2$yhat))))
  print(comp_time)
}