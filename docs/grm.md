# Genomic Relationship Matrix
For a SNP matrix $Z\in \{0,1,2\}^{n \times s}$, the classical genomic relationship matrix (GRM) with uniform scaling introduced by VanRaden 2008 is defined by
$$ G= \frac{PZ Z^TP^T}{2 p^T(1-p)}, \quad \text{with} \quad P = I - \frac{1}{n}1_n 1_n^T, p = \frac{1}{2n}1_n^TZ $$

As noted in Schlather 2020, the computation of $G$ is vastly more efficient when decomposing it into

$$ n^2 2 p^T(1-p) G =n^2 M - n 1_n1_n^TM - nM 1_n 1_n^T + 1_n1_n^TM 1_n 1_n^T, $$
for $M= ZZ^T$, as the expression on the RHS is integer-valued and can be obtained from $M$ with matrix-vector multiplications and rank-k updates. Denoting the columnwise sum of $M$ by $m = M1_n$ and using the symmetry of $M$, the above expression can be reformulated into  
$$  n^2 2 p^T(1-p) G = n^2 M + (-n) ( 1_n m^T + m1_n^T) + (1_n^Tm) 1_n 1_n^T,$$
which corresponds to a `syr2` operation followed by a `syr`.

In practice, $n^2M$ could leave the `Int32` range when $n^2 s > 2^{29}\approx 536$ mn. For the lack of appropriate accumulators on GPUs, we hence convert the entries of $M$ to double-precision floating point values after computation of $M$, as there are hardly any integer-valued BLAS Level-2 operations anyway. 


### References