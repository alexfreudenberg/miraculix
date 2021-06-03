## Authors 
## Martin Schlather, schlather@math.uni-mannheim.de
##
##
## Copyright (C) 2017 -- 2019 Martin Schlather
#               2014 Florian Skene
##
## This program is free software; you can redistribute it and/or
## modify it under the terms of the GNU General Public License
## as published by the Free Software Foundation; either version 3
## of the License, or (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; if not, write to the Free Software
## Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.  


print_haplogeno <- function(x, name, ...) {
  info <- attr(x, "information")
#  print.default(x)
#  print.default(info)
  indiv <- info[INDIVIDUALS + 1]
  method <- info[METHOD + 1]
  if (indiv > 1) {
    cat(if (method == Haplo) "Loci" else "SNP",
        " x Indiv ", name, " matrix with ",
        SNPCODING_NAMES[method + 1],
        " snpcoding: [1:", info[SNPS + 1], ", 1:", indiv,"]",
        sep="")
    if (method == Haplo) cat(" (twice)")
  } else
    cat(if (method == Haplo) "Haplotype " else "SNP",
        "information at ", info[SNPS + 1], " loci.", sep="")
    
  if (is.character(x))
    cat("\n<", x, "> (",
        if (length(attr(x, "M")) == 0) "unloaded" else "loaded", ")", sep="")
  else if (length(f <- attr(x, "filename")) > 0)
    cat("\n<", f, "> (loaded)", sep="")
  cat("\n")
}


as.matrix.haplomatrix <- function(x, ...) decodeHaplo(x, ...)


print.haplomatrix <- function(x, ...)  print_haplogeno(x, "haplo", ...)


str.haplomatrix <- function(object, ...) {
  print.haplomatrix(object, ...)
  cat("Attributes are a ")
  utils::str(base::attributes(object)) ##
}


as.matrix.genomicmatrix <- function(x, ...) decodeGeno(x, ...)
  
print.genomicmatrix <- function(x, ...) print_haplogeno(x, "genomic", ...)

str.genomicmatrix <- function(object, ...) {
  print.genomicmatrix(object, ...)
  cat("Attributes are a ")
  utils::str(base::attributes(object)) ##
  ##class(object) <- NULL;    str(object)
} 
  

as.haplomatrix <- function(object, ...) haplomatrix(object, ...)

as.genomicmatrix <- function(object, ...) genomicmatrix(object, ..., loading = FALSE)
#"as.haplomatrix<-" <- function(object, ...) haplomatrix(object, ..., value)
#  do.call(paste0("as.", value), list(object, ...))
#"as.haplomatrix<-" <- function(object, ...) genomicmatrix(object, ..., value)
#  do.call(paste0("as.", value), list(object, ...))


