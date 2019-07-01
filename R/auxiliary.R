
## Authors 
## Martin Schlather, schlather@math.uni-mannheim.de
##
##
## Copyright (C) 2019 -- 2019 Martin Schlather
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

has.instruction.set <- function(which=c("SSE2", "SSSE3",  "AVX", "AVX2")) {
  ans <- logical(length(which))
  for (i in 1:length(ans)) {
    ans[i] <-
      switch(which[i],
             "SSE2" = .Call(C_hasSSE2),
             "SSSE3" = .Call(C_hasSSSE3),
             "AVX" = .Call(C_hasAVX),
             "AVX2" = .Call(C_hasAVX2),
             FALSE)
  }
  ans
}



#########################################################
## from below here are only temporary functions

unlock <- function(...,  envir = as.environment(-1)) {
  L <- list(...)
  if (!all(sapply(L, function(x) length(attr(x, "information")) > 0)))
    stop("not all arguments are coded objects (Note that 'codeOrigin' formally does return a coded object. So these objects cannot be unlock-ed.");
 if (is.character(L[[1]])) {
    if (!all(sapply(L, is.character)))
      stop("arguments must be all variables or all character vectors")
    L <- unlist(L)
  } else {
    L <- as.character(substitute(list(...)))[-1]
  }
  for (i in L)  .Call(C_unlock, get(i, envir=envir))
}


dolocking <-function(do=TRUE) {
  if (length(do) > 0 && do) stop("sorry dolocking has not been programmed yet")
  .Call(C_dolocking, do)
}


## only for development of rmCoded needed
## check: how to remove an object that has been defined on the upper level
## no further use
RM <- function(...) {
  env <- parent.env(environment())
  L <- list(...)
  if (is.character(L[[1]])) {
    if (!all(sapply(L, is.character)))
      stop("arguments must be all variables or all character vectors")
    L <- unlist(L)
  } else {
    L <- as.character(substitute(list(...)))[-1]
  }
  rm(list=L, envir=env)
}

SAVE <- function(name, ans) {
  print(name) ##
  variable.name <- "__Z"
  ENV <- .GlobalEnv
  name <- paste(variable.name, name, sep=".")
  name <- name[length(name)]
  if (!exists(name, envir=ENV)) assign(name, 0, envir=ENV)
  cur <- get(name, envir=ENV) + 1
  file <- paste(name, cur, "rda", sep=".")
  Print(file, ans) ##
  assign(name, cur, envir=ENV)
  save (file=file, ans) 
  return(ans)
}

RETURN <- function(name) {
  variable.name <- "__Z"
  ENV <- .GlobalEnv
  name <- paste(variable.name, name, sep=".")
  name <- name[length(name)]
  if (!exists(name, envir=ENV)) assign(name, 0, envir=ENV)
  cur <- get(name, envir=ENV) + 1
  file <- paste(name, cur, "rda", sep=".")  
  if (file.exists(file)) {
    cat("playback", file, "\n")
    load (file) 
    assign(name, cur, envir=ENV)
    return(name)
  }
  stop("hierher und nicht weiter")
}

PB <- function(name) {
  variable.name <- "__Z"
  ENV <- .GlobalEnv
  name <- paste(variable.name, name, sep=".")
  name <- name[length(name)]
  if (!exists(name, envir=ENV)) assign(name, 0, envir=ENV)
  cur <- get(name, envir=ENV) + 1
  file <- paste(name, cur, "rda", sep=".")  
  return(file.exists(file))
}


DEBUG <- function(name, population, gen, sex, nr) {
  variable.name <- "__Z"
  ENV <- .GlobalEnv
  name <- paste(variable.name, name, sep=".")
  name <- name[length(name)]
  if (!exists(name, envir=ENV)) assign(name, 0, envir=ENV)
  cur <- get(name, envir=ENV) + 1
  file <- paste(name, cur, "rda", sep=".")  
  if (cur == 2866) {
    .Call(C_Debug);
##    str(population$breeding[[1]][[1]][[318]])
##    str(population$breeding[[gen]][[sex]][[nr]])
##    print(population$breeding[[gen]][[sex]][[nr]][1:2])
  } # else .Call(C_StopDebug);
}

