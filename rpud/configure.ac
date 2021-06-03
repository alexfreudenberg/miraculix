
AC_INIT([miraculix], 0.9)
CXX=`"${R_HOME}/bin/R" CMD config CXX`
AC_PROG_CXX
AC_LANG(C++)

AC_MSG_CHECKING([which configure args were set by the user])

##
if test "${CXX_FLAGS}" ; then
  AC_MSG_RESULT([$CXX_FLAGS])
else
  AC_MSG_RESULT([none])
fi



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
  CPPFLAGS="-m$SET"
  info=`echo "${INFO_CPU}" | grep " $SET "`
  if test "${info}" != ""; then
     AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[]])], [ FX="$FX $CPPFLAGS" ], [ ])
  fi
done

CPPFLAGS="$FX"

##CPPFLAGS="" ## to be deleted if the above does not cause any problems
if test "$CPPFLAGS" != "" ; then
  AC_MSG_RESULT([$CPPFLAGS])
else
  AC_MSG_RESULT([none])
fi


if test "${USE_GPU}" ; then
      CXX_FLAGS="$CXX_FLAGS -DUSEGPU"
      AC_MSG_CHECKING([whether CUDA_HOME is set])
      if test -z "${CUDA_HOME}"; then
          AC_MSG_RESULT(CUDA_HOME not set; using /usr/local/cuda)
          CUDA_HOME="/usr/local/cuda"
      else
          AC_MSG_RESULT(using CUDA_HOME=${CUDA_HOME})
      fi
      AC_CHECK_FILE([${CUDA_HOME}/bin/nvcc], CUDA_INSTALLED=true,)
      if test -z "${CUDA_INSTALLED}"; then
          AC_MSG_RESULT([No CUDA installation found, install CUDA or specify CUDA_HOME.])
      else
          AC_MSG_CHECKING(whether this is the 64 bit linux version of CUDA)
          AC_CHECK_FILE([${CUDA_HOME}/lib64/libcublas.so], , 
          AC_MSG_ERROR([this package only works with 64 bit installations of CUDA]))
          AC_SUBST(CUDA_INSTALLED)

      fi
fi
CPICFLAGS=`"${R_HOME}/bin/R" CMD config CPICFLAGS`
if test "${R_LIB}"; then
    R_LIB="${R_LIB}"
else
    R_LIB="{R_HOME}/lib"
fi

if test -z "${R_INCLUDE_DIR}"; then
    R_INCLUDE_DIR=`"${R_HOME}/bin/R" CMD config --cppflags`
else
    R_INCLUDE_DIR="-I${R_INCLUDE_DIR}"
fi
AC_MSG_RESULT(R_LIBS_USER ${R_LIBS_USER})

AC_MSG_RESULT(using ${R_INCLUDE_DIR} for R header files)
AC_MSG_RESULT(using ${R_LIB} for R shared libraries)

MIRACULIX_CXXFLAGS="$CXX_FLAGS $CPPFLAGS"
AC_SUBST(MIRACULIX_CXXFLAGS)
AC_SUBST(CUDA_HOME)
AC_SUBST(R_LIB)

AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT
