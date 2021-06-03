## This script benchmarks RFutils Cholesky decomposition against R.
library(rgpu)
library(RandomFieldsUtils)

# Matrix sizes
dimensions <- seq(5e3, 30e3, length.out=10)

for(n in dimensions){
    #Simulate positive definite matrix and random vector to solve for
    x <- 1:n
    y <- runif(n)
    M <- exp(-as.matrix(dist(x) / n)) 
    b0 <- matrix(nr=n, runif(n * 1))
    b <- as.matrix(M %*% b0 + runif(n))

    # Warm-up run
    {solvex(M[1:100,1:100]);NULL}
    # Benchmark miraculix
    print(system.time(z <- solvex(M, b)))
    print(range(b - M %*% z))

    # Benchmark base R
    print(system.time(z <- solve(M, b)))
    print(range(b - M %*% z))
}