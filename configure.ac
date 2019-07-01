
AC_INIT([miraculix], 0.6)

USE_FLAGS=true

if "${R_HOME}/bin/R" CMD config CXX11 > /dev/null ; then
  MIRACULIX_CXX=`"${R_HOME}/bin/R" CMD config CXX11`
elif "${R_HOME}/bin/R" CMD config CXX1X > /dev/null ; then
  MIRACULIX_CXX=`"${R_HOME}/bin/R" CMD config CXX1X`
else
  MIRACULIX_CXX=`"${R_HOME}/bin/R" CMD config CXX`
fi

if test "x$TARGET" = "x" ; then
  case `basename $MIRACULIX_CXX` in
  g++ | clang++ )
    if uname -m | grep  -iq sparc ; then
      echo Since this is a sparc processor we set TARGET=default
      TARGET=default
    else
      TARGET=native
    fi
    ;;
  *)
    echo unknown compiler: $MIRACULIX_CXX therefore we set TARGET=default
    TARGET=default
    ;;
  esac
fi

if test "x$MIRACULIX_TARGET" != "x" ; then
  TARGET=$MIRACULIX_TARGET
fi

if test $USE_FLAGS = false ; then
  TARGET=default
fi

echo Using target: $TARGET 
echo Using further args: $@


if test $TARGET = native ; then
  CXX_FLAGS="-O2 -march=native $CXX_FLAGS"
elif test $TARGET = default ; then
  CXX_FLAGS="$CXX_FLAGS"
elif test $TARGET = empty ; then
  echo No compiler flags used
fi


echo Config in `pwd`
echo Config: using CXX_FLAGS=$CXX_FLAGS

MIRACULIX_CXXFLAGS=$CXX_FLAGS

AC_SUBST(MIRACULIX_CXXFLAGS)
AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT
