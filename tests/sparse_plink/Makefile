
CUDA=false
COMPILER=intel

ifeq ($(COMPILER), gcc)
 FORTRAN=gfortran
 FFLAGS=-O3 -ftree-vectorize -mtune=native -march=native 
 FFLAGS+=-cpp -fopenmp 
 LDFLAGS=-llapack -lblas -lgfortran -lm -lquadmath -fopenmp

else ifeq ($(COMPILER), intel)
 FORTRAN=ifort
 FFLAGS=-O3 -mavx2 -mtune=native -march=native
 FFLAGS+=-I${MKLROOT}/include -I${MKLROOT}/include/intel64/lp64 -L${MKLROOT}/lib/intel64
 FFLAGS+=-fpp -qopenmp -qopt-matmul
 LDFLAGS=-liomp5 -lmkl_blas95_lp64 -lmkl_lapack95_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core  -lpthread -lm -ldl -lstdc++
 LDFLAGS+=-g -check all -traceback

endif

SRC_DIR=$(PWD)/../../src
FFLAGS+=-I$(SRC_DIR)/bindings/Fortran
LDFLAGS+=-L$(SRC_DIR)/miraculix -L$(SRC_DIR)/bindings/Fortran -L$(SRC_DIR)/cuda
LDFLAGS+= -Wl,-rpath=$(SRC_DIR)/miraculix -Wl,-rpath=$(SRC_DIR)/bindings/Fortran -Wl,-rpath=$(SRC_DIR)/cuda 
LDFLAGS+=-lmiraculix -lfortran 


ifeq ($(CUDA), true)
LDFLAGS+=-L${CUDA_HOME}/lib64 -lcudart -lcublas -lcusolver -lcusparse -lcusolver -ldgemm_compressed_cuda
FFLAGS+=-DCUDA
endif


OBJ = $(SRC:.f90=.o)
MODS = $(OBJ:.o=.mod)

LIBFORTRAN = $(SRC_DIR)/bindings/Fortran/libfortran.a
LIBMIRACULIX = $(SRC_DIR)/miraculix/libmiraculix.a

EXE = test_sparse_plink
extension = .out



all: $(EXE)


$(EXE): $(EXE).o $(LIBMIRACULIX) $(LIBFORTRAN)
	$(FORTRAN) $(FFLAGS) $(EXE).o -o $@$(extension) $(LDFLAGS) 

%.o: %.f90
	$(FORTRAN) $(FFLAGS) -c $<

$(EXE).o: $(LIBFORTRAN)

clean:
	$(RM) $(OBJ) $(MODS) $(LIB) $(EXE) $(ALL) *.optrpt

tests: all
	./$(EXE)$(extension)

