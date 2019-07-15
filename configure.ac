
AC_INIT([miraculix], 0.8)

dnl USE_FLAGS=false


if "${R_HOME}/bin/R" CMD config CXX11 > /dev/null ; then
  MIRACULIX_CXX=`"${R_HOME}/bin/R" CMD config CXX11`
elif "${R_HOME}/bin/R" CMD config CXX1X > /dev/null ; then
  MIRACULIX_CXX=`"${R_HOME}/bin/R" CMD config CXX1X`
else
  MIRACULIX_CXX=`"${R_HOME}/bin/R" CMD config CXX`
fi


case `basename $MIRACULIX_CXX` in
  g++ | clang++ )
if test $USE_FLAGS = true ; then
  CXX_FLAGS="-O2 -march=native $CXX_FLAGS"
else
  CXX_FLAGS="-mavx $CXX_FLAGS"
fi
;;
  *)
echo Unknown compiler: $MIRACULIX_CXX, no compiler flags used
    ;;
esac

echo Config: using CXX_FLAGS=$CXX_FLAGS

MIRACULIX_CXXFLAGS=$CXX_FLAGS

AC_SUBST(MIRACULIX_CXXFLAGS)
AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT

