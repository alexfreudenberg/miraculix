

summary.RFopt <- function(object, ...) {  
  object <- lapply(object, function(z) z[order(names(z))])
  object <- object[c(1, 1 + order(names(object[-1])))]
  class(object) <- "summary.RFopt"
  object
}


print.summary.RFopt <- function(x, ...) {
  str(x, give.attr=FALSE, ...) #
  invisible(x)
}

print.RFopt <- function(x, ...) {
  print.summary.RFopt(summary.RFopt(x, ...)) #
  invisible(x)
}

summary.RFoptElmnt <- function(object, ...) {
  object <- object[order(names(object))]
  class(object) <- "summary.RFoptElmt"
  object
}

print.summary.RFoptElmnt <- function(x, ...) {
  str(x, give.attr=FALSE, ...) #
  invisible(x)
}

print.RFoptElmnt <- function(x, ...) {
  print.summary.RFoptElmnt(summary.RFoptElmnt(x, ...)) #
  invisible(x)
}

detach_packages <- function(pkgs) {
  for (pkg in pkgs) {
    pkg <- paste0("package:", pkg)
   while(pkg %in% search()) detach(pkg, unload = TRUE, character.only=TRUE)
  }
}
libraries <- function(pkgs, control, verbose=FALSE) {
  if (length(control) > 0) {
    idx <- pmatch(names(control), names(as.list(args(library))))
    control <- control[idx[!is.na(idx)]]
  }
  for (pkg in pkgs) do.call("library", c(list(pkg), control))
  if (verbose) message("libraries attached.")
}

S <- function(x) if (length(x) > 1) "s" else ""
ARE <- function(x) if (length(x) > 1) "are" else "is"
HAVE <- function(x) if (length(x) > 1) "have" else "has"



reinstallPackages <- function(ic, basic, install.control) {    
  install <- basic$install
##  Print(install, interactive(), ic)
  verbose <- TRUE
  force <- FALSE
  quiet <- FALSE
  if (ic) {
    if (length(install.control$verbose)) verbose <- install.control$verbose
    if (is.logical(install.control$quiet)) quiet <- install.control$quiet      
    if (install %in% c("ask", "none")) install <- "install"
    if (is.logical(install.control$force)) {
      force <- install.control$force
      if (!force) {
        ##       if (!interactive())
        ##         stop("'force=FALSE' must be run in the interactive mode of R")
        install <- "ask"
      }
    }
  }

  #Print("A", install.control)
  
  pkgs <- .Call(C_getPackagesToBeInstalled, force)
  if (ic && length((install.control$pkgs))) {
    pkgs <-install.control$pkg
    install.control$pkgs <- NULL
  }
  
  verbose <- verbose && !quiet
  if (length(pkgs) == 0) {
    ## should/can happen only if user calls it
    ## Print(basic, basic$printlevel)
    .Call(C_AVXmessages, "all")
    if (verbose)
      message("No packages found to be installed. This happens particularly if the the installation process was interrupted. Try it again in the next session or use 'RFoptions(install.control=list(force=TRUE))' for instance.")
    return()
  }
  if (install == "ask") {
    omp <- .Call(C_AVXmessages, pkgs)
    if (!quiet)
      cat("The package", S(pkgs), " ", paste0("'", pkgs, "'", collapse=", "),
          " ",
          HAVE(pkgs), " been compiled without appropriate SIMD/AVX2 flags. So, calculations can be slow. If the package",S(pkgs), " ", ARE(pkgs), " recompiled with the necessary flags, the calculations might be faster.\nR should be restarted after re-compiling. The argument 'install.control' might be used to run the re-compilation without asking and to pass further arguments to 'install.packages', e.g., 'RFoptions(install.control=list(repos=NULL, verbose=TRUE))'\nTo avoid this feedback, set 'RFoptions(install=\"none\")' or 'RFoptions(install=\"install\")' before calling any other function of '", pkgs[length(pkgs)],"'.\n", sep="")
    txt <- paste0("Shall the package", S(pkgs),
                  " and all further packages based on 'RandomFieldsUtils' be recompiled (Y/n/args) ? ")
    install.control <- readline(txt)
    install <- if (install.control %in% c("n", "N")) "none" else "install"
    if (nchar(install.control) <= 3)  install.control <-""
    if (verbose) {
      if (install == "none") {
        cat("If you have stopped the re-compilation since does not work, consider one of the following possiblities:")
        .Call(C_AVXmessages, NULL)
        cat("\nIf all fails, call 'RFoptions(install=\"none\")' after any loading of the package.\nOtherwise you will be bothered with the above question again and again.\n")
      } else {
        S <- "\t*************************************************\n"
        cat("\n", S, "\t***         Do not forget to restart R.       ***\n",S)
        sleep.milli(1500)
      }
    }
  } else omp <- .Call(C_AVXmessages, "OMP")
  ##      Print(pkgs, install.control, install)
  if (install != "none") {
    if (is.character(install.control)) 
      install.control <- eval(parse(text=paste("list(", install.control, ")")))
    CXX_FLAGS <- args <- ""
    if (length(install.control$configue.args) > 0) {
      args <- install.control$configue.args
      install.control$configue.args <- NULL
    }
    if (length(install.control$CXX_FLAGS) > 0) {
      CXX_FLAGS <- install.control$CXX_FLAGS
      install.control$CXX_FLAGS <- NULL
    }
    idx <- pmatch(names(install.control),names(as.list(args(install.packages))))
    ##        Print(idx, names(install.control), names(install.control))        
    install.control <- install.control[which(!is.na(idx))]
    
    args <- paste0(args, " CXX_FLAGS='", CXX_FLAGS, " ", if (nchar(omp) > 0) omp
                 , "'")
    
    if (verbose) Print(install.control, args) ## OK
    args <- paste0("USE_AVX='yes' TRY_GPU='yes' ", args)
    for (pkg in pkgs) {
      z <- Try(do.call("install.packages",
                       c(list(pkgs=pkg, type="source", configure.args=args),
                         install.control)))
      if (is(z, "try-error")) print(z) ## OK
    }
    ## on.exit({detach_packages(rev(pkgs)); libraries(pkgs)}, add=TRUE)
  }
  cat("\n\n")
}


RFoptions <- function(..., no.readonly=TRUE, no.class=FALSE, install.control=NULL) {
  opt <- .External(C_RFoptions, ...)  
  ##  if (is.list(opt)) Print(basic) else Print(opt)
  ic <- hasArg("install.control")
  if (ic || (length(opt) > 0 && is.list(opt$basic) &&
             opt$basic$installPackages && interactive())) {
    reinstallPackages(ic=ic, basic=opt$basic, install.control=install.control)
    if (ic) return(invisible(NULL))
  }
  if (length(opt) == 0 || no.class) return(invisible(opt))
  if (is.list(opt[[1]])) {
    opt <- lapply(opt,
		  function(x) {
		    class(x) <- "RFoptElmnt"
		    x
		})
    class(opt) <-  "RFopt"
  } else class(opt) <-  "RFoptElmnt"
  if (!no.readonly) {
    opt$readonly <- list()
  }
  opt
}
