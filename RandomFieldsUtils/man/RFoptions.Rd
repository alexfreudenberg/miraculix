\name{RFoptions}
\alias{RFoptions}
\alias{PIVOT_NONE}
\alias{PIVOT_AUTO}
\alias{PIVOT_DO}
\alias{PIVOT_IDX}
\alias{LA_AUTO}
\alias{LA_INTERN}
\alias{LA_R}
\alias{LA_GPU}
%\alias{PIVOT_IDXBACK}
\alias{PIVOTSPARSE_MMD}
\alias{PIVOTSPARSE_RCM}
\title{Setting control arguments}
\description{
 \command{\link{RFoptions}} sets and returns control arguments for the analysis
 and the simulation of random fields
}
\usage{
 RFoptions(..., no.readonly = TRUE, no.class=FALSE, install.control=NULL)
}
\arguments{
  \item{...}{arguments in \code{tag = value} form, or a list of tagged
    values.}
  \item{no.readonly}{If \command{\link{RFoptions}} is called without
    argument then all arguments are returned in a list. If
    \code{no.readonly=TRUE} then only rewritable arguments are returned.
    Currently all arguments are rewritable. So, no effect for package
    \pkg{RandomFieldsUtils}
  }
  \item{no.class}{logical. If \code{TRUE} the list is returned without
    class specification.
  }
  \item{install.control}{list, containing any argument 
    \command{\link[utils]{install.packages}} except \code{type}.
    The argument \code{configure.args} (i) is merged with the presetting
    of \command{RFoptions}; (ii) may not contain \code{'CXX_FLAGS'}
    which should be passed as an extra argument with the list.
    
    Note that if this argument is given (even with value \code{NULL}),
    an immediate installation takes place. In case the user tries to
    force to install 0 packages, an overview over the packages is given.
    If the user is asked whether re-installation shall take place,
    user can pass arguments to install.packages, e.g.,
    \code{"repos=NULL, quiet=TRUE"}.
    
    If this argument is given, no further
    argument may be passed to \command{RFoptions}.
    
    The arguments \code{verbose} and \code{quiet} also affect\code{RFoptions}.
    
    Additionally, the argument \code{force} has the following behaviour:
    \itemize{
      \item{\code{TRUE}}{
	reinstallation of all attached libraries based on and including
	\pkg{RandomFieldsUtils}
      }
      \item{\code{FALSE}}{In case some packages have to be re-installed
	the user will be asked.
      }
      \item{not given}{
	reinstallation of the attached libraries based on and including
	\pkg{RandomFieldsUtils} that have not been tried yet in the
	current session.
      }
    }
  }
}
\details{
  The subsections below comment on\cr
  \bold{1. \code{basic}: Basic options}\cr
  \bold{2. \code{solve}: Options for solving linear systems}\cr
  \bold{3. Reserved words}\cr
  \cr

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % \bold{16. Options for RFloglikelihood}\cr
  % 
  % "auto", "full", "composite", "selection"
  

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  \bold{1. Basic options}
  \describe{
    \item{\code{asList}}{logical. Lists of arguments are treated slightly
      different from non-lists. If \code{asList=FALSE} they are treated the
      same way as non-lists. This options being set to \code{FALSE} after
      calling \command{RFoptions} it should be set as first element of a list.
      
      Default: \code{TRUE}
    }
    
    \item{\code{cores}}{
      Number of cores for multicore algorithms;
      currently only used for the Cholesky decomposition.
      
      Default : \code{1} if the package has been compiled with
      standard flags of CRAN and
      \code{0.75 * cores() + 0.25 * cpus()} else.
    }
    
    \item{\code{cPrintlevel}}{
      \code{cPrintlevel} is automatically set to \code{printlevel}
      when \code{printlevel} is changed.
      Standard users will never use a value higher than 3.
      
      0 : no messages\cr
      1 : messages and warnings when the user's input looks odd\cr
      2 : messages (and internal errors) documenting the choice of the
      simulation method\cr
      3 : further user relevant informations\cr
      4 : information on recursive function calls\cr
      5 : function flow information of central functions \cr
      6 : errors that are internally treated\cr
      7 : details on building up the covariance structure\cr
      8 : details on taking the square root of the covariance matrix\cr
      9 : details on intermediate calculations\cr
      10 : further details on intermediate calculations\cr
      
      Note that \code{printlevel} works
      on the R level whereas \code{cPrintlevel} works on the C level.
      
      Default: 1 \cr % [also  do].\cr
    }

  
    \item{\code{helpinfo}}{logical. If \code{TRUE} then
      additional information is printed for more efficient programming
      in R.

      Default: \code{TRUE}
    }
    
    \item{\code{install}}{character. Only used by linux systems and
      alike including macOS
     The default by CRAN is that SIMD/AVX cannot be used at full extend.
      \code{install} determines what the action if the compiled
      version does not use the full CPU capacities.
      The users usually use 
      \describe{
	\item{\code{"none"}}{no re-installation}
	\item{\code{"ask"}}{asks whether the library should be reinstalled,
	  using the full capacity of the CPU according to the package.}
	\item{\code{"install"}}{install without asking}
      }
      Note that, in each session, a package can be reinstalled only.
      This feature avoids trying to run jobs that cannot be done (e.g.\ due to
      missing programs of the OS). See argument \code{install.control}
      for overwriting this behaviour.
      
      Default: at starting point it is \code{"ask"} or \code{"none"},
      but the value may change during the session.
      
%      The user can also indicate what whould be installed:
%      \describe{
%	\item{\code{"sse"}}{installs SSE}
%	\item{\dots}
%	\item{\code{"avx2"}}{installs AVX2}
%      }

    }

    \item{\code{installPackages}}{logical. Read only.
      Indicates whether packages are left to be re-installed.
      \code{RFoptions(install="none")} sets it to \code{FALSE}.
      \code{RFoptions(install="none", install="ask")} sets it to \code{TRUE}.
    }
    
    \item{\code{kahanCorrection}}{obsolete.
      logical. If \code{TRUE}, the Kahan summation algorithm is used for
      calculating scalar products.
      
      Default: false
    }
   
   \item{\code{la_mode}}{determines
      \describe{
	\item{LA_AUTO, \code{"auto"}}{
	  If a graphic card is recognized, \code{LA_GPU} is used.
	  In all other cases the default is primarily \code{LA_R}.
	  Only on linux systems, the package peforms a simple speed test and
	  takes \code{LA_INTERN} if it is faster than {LA_R}; the time,
	  hence the choice, depends also on the number of cores used.}
  	\item{LA_INTERN, \code{"intern"}}{mostly own algorithms, often
	  based on SIMD/AVX. This option is of
	  interest only if no advanced BLAS/LAPACK
	  has been compiled into R}
  	\item{LA_R, \code{"R"}}{BLAS/LAPACK implementation used by R}
 	\item{LA_GPU, \code{"GPU"}}{
	  This option is available when the package has been compiled with nvcc.
	}
      }
	
      Default: \code{LA_AUTO}
    }

    \item{\code{printlevel}}{If \code{printlevel}\eqn{\le0}{<=0}
      there is not any output on the screen. The
      higher the number the more tracing information is given.
      Standard users will never use a value higher than 3.
      
      0 : no messages\cr
      1 : important (error) messages and warnings\cr
      2 : less important messages\cr
      3 : details, but still for the user\cr
      4 : recursive call tracing\cr
      5 : function flow information of large functions\cr
      6 : errors that are internally treated\cr
      7 : details on intermediate calculations\cr
      8 : further details on intermediate calculations\cr
      
      Default: 1 %[also do].\cr
    }
    
    \item{\code{seed}}{integer (currently only used by the package
      RandomFields). 
      If \code{NULL} or \code{NA}
   %   \command{\link[base]{set.seed}}
       \command{set.seed}
      is \bold{not} called.
      Otherwise,
  %    \code{\link[base]{set.seed}(seed)}
      \code{set.seed(seed)}
      is set
      before any simulations are performed.
      
      If the argument is set locally, i.e., within a function,
      it has the usual local effect. If it is set globally, i.e. by
      \command{RFoptions} the \code{seed} is fixed
      for \bold{all subsequent} calls. 
      
      If the number of simulations \code{n} is greater than one
      and if \code{RFoptions(seed=seed)} is set, the \eqn{i}th
      simulation is started with the seed \sQuote{\code{seed}\eqn{+i-1}}.
      %   The function \code{set.seed} should not be used in case \code{n}
      %    is greater than 1.
      %
      %Vgle!   
      %set.seed(5)
      %RFsimulate(RPschlather(RMmatern(nu=2), xi=1, mu=1, s=1), x, grid=F, n=5)@data
      %set.seed(5)
      %RFsimulate(RPschlather(RMmatern(nu=2.01), xi=1, mu=1, s=1), x,grid=F,n=5)@data
      %RFoptions(cPr=3, seed=5)
      %RFsimulate(RPschlather(RMmatern(nu=2), xi=1, mu=1, s=1), x, grid=F, n=5)@data
      %RFsimulate(RPschlather(RMmatern(nu=2.01), xi=1, mu=1, s=1), x, grid=F,n=5)@data
      
    }
    
    \item{\code{skipchecks}}{logical.
      If \code{TRUE}, several checks whether the given parameter values
      and the dimension are within the allowed range is skipped.
      Do not change the value of this variable except you really
      know what you do.
      
      Default: \code{FALSE}  %[also do]. 
    }
    
    \item{\code{verbose}}{logical. If \code{FALSE} it identical to
      \code{printlevel = 1} else to \code{printlevel = 2}.
    }


    \item{\code{warn_parallel}}{Logical.
      \pkg{RandomFieldsUtils} and
      packages using it, such as \pkg{RandomFields} and \pkg{miraculix},
      should now be prepared for parallelization using package
      \code{parallel}, for instance. Internal OMP parallelization
      of \pkg{RandomFieldsUtils} is done, but only at a view points
      of the subsequent packages.

      As a few parts cannot be in parallel technically
      or from a logical point of view, a hint or a
      warning is given, if such a point is not accessed adequately.
      These messages can be turned off by \code{warn_parallel = FALSE}.      
      
      Default: \code{TRUE}.
    }

  \item{\code{warn_unknown_option}}{integer.       
      \describe{
	\item{\code{0},\code{1},\code{-1}}{
	  Unknown options are all ignored. If the value is positive,
	  a hint is delivered whenever an unknown option is ignored.
	}
	\item{\code{-2},\code{2}}{
	  Unknown options that start with a capital letter are ignored.
	  All others lead to an error. (Note that all \code{RFoptions}
	  start with a minor letter.)
	  If the value is positive,
	  a hint is delivered whenever an unknown option is ignored.
	}
	\item{\code{3},\code{-3}}{
	  Unknown options that consists of a single capital letter are ignored.
	  All others lead to an error. (Note that all \code{RFoptions}
	  start with a minor letter.)
	  If the value is positive,
	  a hint is delivered whenever an unknown option is ignored.	  
	}
	\item{\code{4}}{(and other values) Any unknown option leads to
	  an error.}
      }

      Default for \pkg{RandomFieldsUtils}: \code{3}

      Default for \pkg{RandomFields}: \code{1}
      
    }
  }

    
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  \bold{2. \code{solve}: Options for solving linear systems}
  
  
  \describe{ % to do
    \item{\code{det_as_log}}{
    }
    \item{\code{eigen2zero}}{
      When the svd or eigen decomposition is calculated,
      all values with modulus less than or equal to \code{eigen2zero}
      are set to zero.
      
      Default: \code{1e-12}
    }

    \item{\code{max_chol}}{integer. Maximum number of rows of a matrix in
      a Cholesky decomposition
      
      Default: \eqn{16384}
    }
    \item{\code{max_svd}}{integer. Maximum number of rows of a matrix in
      a svd decomposition
      
      Default: \eqn{10000}
    }

    \item{\code{pivot_partialdet}}{logical. If \code{TRUE} then in case
      of low-rank matrices the determinant is calculated only in the part
      with positive eigenvalues}
      
    \item{\code{pivot}}{Type of pivoting for the Cholesky
      decomposition. Possible values are 
      \describe{
	\item{PIVOT_NONE, \code{"none"}}{No pivoting.}
	\item{PIVOT_AUTO, \code{"auto"}}{If the  matrix has a size greater than
	  3x3 and Choleskey fails without pivoting, privoting
	  is done. For matrices of size less than 4x4, no pivoting and
	  no checks are performed. See also \code{PIVOT_DO}}
 	\item{PIVOT_DO, \code{"do"}}{Do always pivoting.
	  NOTE: privoted Cholesky decomposition yields only very approximately
	  an upper triangular matrix L, but still L^t L = M holds true.
%	  The information about the pivoting sequence are attributed to
%	  the result.
	}
 	\item{PIVOT_IDX, \code{"idx"}}{uses the same pivoting as in the previous
	  pivoted decomposition. This option becomes relevant only when
	  simulations with different parameters or different models shall be
	  performed with the same seed so that also the pivoting must be
	  coupled.
%	  The information about the pivoting sequence are attributed to
%	  the result.
	}
 %	\item{PIVOT_IDXBACK}{ same as \code{PIVOT_IDX}, but
%	the sequence of indices of the pivoting is returned via
%	\code{RFoptions()$solve$pivot_idx}.}
      }

      Default: \code{PIVOT_NONE} 
      
    }
    
    \item{\code{pivot_actual_size}}{integer.
      Genuine dimension of the linear mapping given by a matrix in \link{cholx}.
      This is a very rarely used option when pivoting with
      \code{pivot=PIVOT_IDX}.
    }
     
    \item{\code{pivot_check}}{logical. Only used in pivoted Cholesky
      decomposition. 
      If \code{TRUE} and a numerically zero diagonal element is detected,
      it is checked whether the offdiagonal elements are numerically zero 
      as well.
      (See also \code{pivot_max_deviation} and
      \code{pivot_max_reldeviation}.) 
      If \code{NA} then only a warning is given.
      
      Default: \code{TRUE}
   }

       \item{\code{pivot_idx}}{vector of integer.
      Sequence of pivoting indices in pivoted Cholesky decomposition.
      Note that
      \code{pivot_idx[1]} gives the number of indices that will be
      used. The vector must have at least the length
      \code{pivot_idx[1] + 1}.

      Default: \code{NULL}
    }
    
    \item{\code{pivot_relerror}}{positive number.
      Tolerance for (numerically) negative eigenvalues and for (numerically)
      overdetermined systems appearing in the pivoted Cholesky decomposition.

      Default: \code{1e-11}
    }
    
    \item{\code{pivot_max_deviation}}{positive number.
      Together with \code{pivot_max_reldeviation} it determines
      when the rest of the matrix (eigenvalues) in the pivoted
      Cholesky decomposition are considered as zero.
      
      Default: \code{1e-10}
    }

    \item{\code{pseudoinverse}}{logical. In case of a singular matrix \eqn{M},
      shall the pseudo inverse be returned for \code{solvex(M)}?
      
      Default: \code{FALSE}
    }
    
    \item{\code{pivot_max_reldeviation}}{positive number.
	Together with \code{pivot_max_deviation} it determines
      when the rest of the matrix (eigenvalues) in the pivoted
      Cholesky decomposition are considered as zero.
 
      Default: \code{1e-10}
   }

     
   \item{\code{solve_method}}{
      vector of at most 3 integers that gives the sequence of methods
      in order to inverse a matrix or to calculate its square root:
      \code{"cholesky"}, \code{"svd"}, \code{"lu"}, % \code{"qr"},
      \code{"eigen"} \code{"sparse"},
      \code{"method undefined"}. In the latter case, the algorithm decides
      which method might suit best.

    %  Note that \code{"qr"} returns only the transposed of the inverse!
      
      Note that if \code{use_spam} is not \code{false} the algorithm
      checks whether a sparse matrix algorithm should be used and which is
      then tried first.
      
      Default: \code{"method undefined"}.
    }
    \item{\code{spam_factor}}{
      integer. See argument \code{spam_sample_n}.
      
      Default: 4294967
    }
    
    \item{\code{spam_min_n}}{
      integer vector of size 2. The minimal size for a matrix to apply a
      sparse matrix algorithms automatically. The second value is used
      in case the GPU is activated.
      
      Default: \code{c(400, 4000)}
    }
    
    \item{\code{spam_min_p}}{
      number in \eqn{(0,1)} giving the proportion of
      zero above which an sparse matrix algorithm is used.
      The second value is used
      in case the GPU is activated.
      
      Default: \code{c(0.8, 0.8)}
    }
    \item{\code{spam_pivot}}{
      integer. Pivoting algorithm for sparse matrices:
      \describe{
	\item{PIVOT_NONE}{No pivoting}
	\item{PIVOTSPARSE_MMD}{}
	\item{PIVOTSPARSE_RCM}{}
      }
      See package \code{spam} for details.
      
      Default: PIVOTSPARSE_MMD
    }
    \item{\code{spam_sample_n}}{
      Whether a matrix is sparse or not is tested by a
      \sQuote{random} sample of size \code{spam_sample_n};
      The selection of the sample is iteratively
      obtained by multiplying the index by \code{spam_factor}
      modulo the size of the matrix.
      
      Default: 500. 
    }
    \item{\code{spam_tol}}{
      largest absolute value being considered as zero.
      Default: \code{DBL_EPSILON}
    }
    \item{\code{svdtol}}{Internal.
      When the svd decomposition is used for calculating the square root of
      a matrix then the absolute componentwise difference between
      this matrix and the square of the square root must be less
      than \code{svdtol}. No check is performed if
      \code{svdtol} is not positive.
      
      Default: \code{0}
    }
     \item{\code{use_spam}}{
      Should the package \code{spam} (sparse matrices)
      be used for matrix calculations?
      If \code{TRUE} \pkg{spam} is always used. If \code{FALSE},
      it is never used. If \code{NA} its use is determined by
      the size and the sparsity of the matrix.
      
      Default: \code{NA}. 
    }
  }

  \bold{3. Reserved Words}
  \describe{
    \item{\code{list_}}{
      \code{list_} usually equals the output of \code{RFoptions()}.
      This argument is used to reset the RFoptions.
      Some of the options behave differently if passed through
      \code{list_}. E.g. a warning counter is not reset.
      The argument \code{list_} cannot be combined with any other arguments.
    }
    \item{\code{getoptions_}}{string vector of prefixes that indicate
      classes of options. In this package they
      can be \code{"basic"} and \code{"solve"}. (E.g. package
      \pkg{RandomFields} has many more classes of options.)
      The given classes of options  are then
      returned by \code{RFoptions()}. Note that the values are the
      previous values.
      
      \code{getoptions_} must always be the very first argument.
    }
    \item{\code{saveoptions_}}{string vector of prefixes. Same as for
      \code{getoptions_}, except that important classes are always
      returned and thus  should not be given. Hence \code{saveoptions_}
      is often a convenient short cut for \code{getoptions_}.
      The class always included in this package is \code{"basic"}, in
      package \pkg{RandomFields} these are the two classes
      \code{"basic"} and \code{"general"}.
   
      \code{saveoptions_} must always be the very first argument. In
      particular, it may not given at the same time with \code{getoptions_}.
    }
    \item{\code{local_}}{logical. This options is allowed only when
      advanced packages are used, see \pkg{RandomFields}.
    }
    \item{\code{warnUnknown_}}{integer. Same as option
      \code{warn_unknown_option}, except that its value overwrites
      the value of \code{warn_unknown_option} in the current command
      \code{RFoptions}. This options must be placed between \code{CODE}
      and \code{getoptions_}, if the latter are used.
    }
  }
}


 
\value{
 \code{NULL} if any argument is given, and the full list of
 arguments, otherwise.
}

\me

\examples{
%   library(RandomFieldsUtils)
n <- 10
M <- matrix(1, ncol=n, nrow=n)

\dontrun{
try(chol(M)) ## error, as M is not strictly positive definite
try(cholx(M)) ## also fails
}

%cholx(M) ## works

RFoptions(la_mode=LA_INTERN, pivot=PIVOT_AUTO)
cholx(M) ## works
RFoptions(la_mode=LA_R)


RFoptions(solve_method="svd", pseudoinverse=TRUE)
solvex(M)
RFoptions(solve_method="method undefined", pseudoinverse=FALSE)

}

\keyword{spatial}

