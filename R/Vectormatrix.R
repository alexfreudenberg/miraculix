


vector012matrix <- function(v, M) .Call(C_vector012matrix, v, M)
matrixvector012 <- function(M, v) .Call(C_matrixvector012, M, v)


## not used yet:
"%*%" <- function(x, y) {
  if (is(x, "genomicmatrix")) {
    stopifnot(is.vector(y))
    matrixvector012(x, y)
  } else if (is(y, "genomicmatrix")) {
    stopifnot(is.vector(x))
    vector012matrix(x, y)
  } else {
    base::"%*%"(x, y)
  }
}

crossprod <- function(x, y=NULL) {
  if  ((xis <- is(x, "genomicmatrix")) || (yis <- is(y,  "genomicmatrix"))) {
    stopifnot(is.null(y))
    crossprodx(SNPxIndiv=x)
  } else {
    base::crossprod(x, y)
  }
}
