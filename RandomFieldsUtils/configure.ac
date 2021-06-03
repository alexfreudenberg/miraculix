
AC_INIT([RandomFieldsUtils], 0.6)
CXX=`"${R_HOME}/bin/R" CMD config CXX`
AC_PROG_CXX
AC_LANG(C++)


AC_MSG_NOTICE([value of USE_AVX is '${USE_AVX}'])
AC_MSG_NOTICE([value of USE_GPU is '${USE_GPU}'])
AUTOCPPFLAGS=""
MY_OWN_CU_FILES=""
CUDA_LINKS=""
if test "${USE_AVX}" ; then
  CPU_FLAGS=""
  AC_CHECK_PROG(LSCPU, [lscpu], yes, no)

  if test "x${LSCPU}" == xyes; then
    INFO_CPU=`lscpu | grep Flags | tr "[[:upper:]]" "[[:lower:]]"`
  else
    AC_CHECK_PROG(SYSCTL, [sysctl], yes, no)
    if test  "x${SYSCTL}" == xyes; then
      INFO_CPU=`sysctl -a 2> /dev/null | grep machdep.cpu.features | tr "[[:upper:]]" "[[:lower:]]"`
    fi
  fi


  AC_MSG_CHECKING([which SIMD flags can be recognized])
  FX=""
  ## sysctl -a | grep "cpu.features:" # for OS X
  for SET in sse2 ssse3 avx2; do 
    AUTOCPPFLAGS="-m$SET"
    info=`echo "${INFO_CPU}" | grep " $SET "`
    if test "${info}" != ""; then
     AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[]])], [ FX="$FX $AUTOCPPFLAGS" ], [ ])
    fi 
  done

  AUTOCPPFLAGS="$FX"   
  if test "$AUTOCPPFLAGS" != "" ; then
    AC_MSG_RESULT([$AUTOCPPFLAGS])
  else	
    AC_MSG_RESULT([none])
  fi
fi


ORIG_CXX_FLAGS="$CXX_FLAGS"

AC_MSG_NOTICE([values of USE_GPU are +${USE_GPU}+${TRY_GPU}+])
if test "${USE_GPU}${TRY_GPU}" ; then
      CXX_FLAGS="$CXX_FLAGS -DUSEGPU"
      AC_MSG_CHECKING([whether CUDA_HOME is set])
      if test -z "${CUDA_HOME}"; then
          AC_MSG_RESULT([CUDA_HOME not set; using /usr/local/cuda])
          CUDA_HOME="/usr/local/cuda"
      else
          AC_MSG_RESULT(using CUDA_HOME=${CUDA_HOME})
      fi
      AC_CHECK_FILE([${CUDA_HOME}/bin/nvcc], CUDA_INSTALLED=true,)
      if test "${CUDA_INSTALLED}"; then
          AC_MSG_CHECKING(whether this is the 64 bit linux version of CUDA)
          AC_CHECK_FILE([${CUDA_HOME}/lib64/libcublas.so], , 
          AC_MSG_ERROR([this package only works with 64 bit installations of CUDA]))
          AC_SUBST(CUDA_INSTALLED)
          MY_OWN_CU_FILES="solve_gpu.o"
	  CUDA_LINKS=" -L${CUDA_HOME}/lib64 -lcudart -lcublas -lcusolver -lcusolverMg"
       else
        if test "${USE_GPU}" ; then
          AC_MSG_ERROR([No CUDA installation found, install CUDA or specify CUDA_HOME.])
        fi
 	CXX_FLAGS="$ORIG_CXX_FLAGS"
      fi
fi



## install.packages("miraculix_0.9.8.tar.gz", configure.args="CXX_FLAGS=-march=native", repos=NULL)
MY_OWN_CXXFLAGS="$CXX_FLAGS $AUTOCPPFLAGS"
MY_OWN_LIBFLAGS="$LIB_FLAGS"

AC_SUBST(MY_OWN_CXXFLAGS)
AC_SUBST(MY_OWN_LIBFLAGS)
AC_SUBST(MY_OWN_CU_FILES)
AC_SUBST(CUDA_HOME)
AC_SUBST(CUDA_LINKS)
AC_SUBST(R_LIB)


AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT
