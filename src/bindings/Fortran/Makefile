COMPILER=intel

ifeq ($(COMPILER), gcc)
 FORTRAN=gfortran
 FFLAGS=-O3 -ftree-vectorize -mtune=native -march=native 
 FFLAGS+=-cpp -fopenmp 
 LINK_FLAGS=-llapack -lblas -lgfortran -lm -lquadmath -fopenmp

 LINKER=-std=gnu++14 -L/usr/local/lib -L/usr/lib64

else ifeq ($(COMPILER), intel)
 FORTRAN=ifort
 FFLAGS=-O3 -mavx2 -mtune=native -march=native
 FFLAGS+=-I${MKLROOT}/include -I${MKLROOT}/include/intel64/lp64 -L${MKLROOT}/lib/intel64
 FFLAGS+=-fpp -qopenmp -qopt-matmul -qopt-report=5
 LINK_FLAGS=-liomp5 -lmkl_blas95_lp64 -lmkl_lapack95_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core  -lpthread -lm -ldl -lstdc++

 LINKER=-std=gnu++14

endif



LINK=$(FORTRAN)


SRC = \
      mod5codesapi.f90 \
      modplink_miraculix.f90 \
      modtestplink.f90

OBJ = $(SRC:.f90=.o)
MODS = $(OBJ:.o=.mod)

LIB = libfortran.a

lib: $(LIB)

$(LIB): $(OBJ)
	ar cr $@ $(OBJ)

%.o: %.f90
	$(FORTRAN) $(FFLAGS) -c $<

clean:
	$(RM) $(OBJ) $(MODS) $(LIB) $(ALL) *.optrpt

.PHONY: clean lib
