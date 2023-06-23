# Genomic Relationship Matrix
For a SNP matrix $Z\in \{0,1,2\}^{n \times s}$, the classical genomic relationship matrix (GRM) with uniform scaling introduced by VanRaden 2008[^1] is defined by
$$G= \frac{PZ Z^TP^T}{2 p^T(1-p)}, \quad \text{with} \quad P = I - \frac{1}{n}1_n 1_n^T, p = \frac{1}{2n}1_n^TZ$$

As noted in Schlather 2020[^2], the computation of $G$ is vastly more efficient when decomposing it into

$$n^2 2 p^T(1-p) G =n^2 M - n 1_n1_n^TM - nM 1_n 1_n^T + 1_n1_n^TM 1_n 1_n^T,$$
for $M= ZZ^T$, as the expression on the RHS is integer-valued and can be obtained from $M$ with matrix-vector multiplications and rank-k updates. Denoting the column-wise sum of $M$ by $m = M1_n$ and using the symmetry of $M$, the above expression can be reformulated into  
$$n^2 2 p^T(1-p) G = n^2 M + (-n) ( 1_n m^T + m1_n^T) + (1_n^Tm) 1_n 1_n^T,$$
which corresponds to a `syr2` operation followed by a `syr`.

In practice, $n^2M$ could leave the `Int32` range when $n^2 s > 2^{29}\approx 536$ mn. For the lack of appropriate accumulators on GPUs, we hence convert the entries of $M$ to double-precision floating-point values after computation of $M$, as there are hardly any integer-valued BLAS Level-2 operations anyway. 

Unfortunately, this approach fails when the standardized GRM $\tilde{G}$ is required, which scales the columns of $PZ$ by standard deviation of each SNP:
$$\tilde{G) = (PZB)(PZB)^T \quad \text{with} \quad B = \text{diag}( 2 p^T (1-p))^{1/2},$$


### References
[^1]: VanRaden PM. Efficient methods to compute genomic predictions. J Dairy Sci. 2008 Nov;91(11):4414-23. doi: 10.3168/jds.2007-0980. PMID: 18946147.

[^2]: Schlather, Martin. Efficient calculation of the genomic relationship matrix. bioRxiv (2020): 2020-01.

[^3]: Mäntysaari EA, Evans RD, Strandén I. Efficient single-step genomic evaluation for a multibreed beef cattle population having many genotyped animals. J Anim Sci. 2017 Nov;95(11):4728-4737. doi: 10.2527/jas2017.1912. PMID: 29293736; PMCID: PMC6292282.
