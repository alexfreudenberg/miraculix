CUDA=true
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

endif

# Paths and references
SRC_DIR=$(PWD)/../../src
FFLAGS+=-I$(SRC_DIR)/bindings/Fortran
LDFLAGS+=-L$(SRC_DIR)/miraculix -L$(SRC_DIR)/bindings/Fortran -L$(SRC_DIR)/cuda
LDFLAGS+= -Wl,-rpath=$(SRC_DIR)/miraculix -Wl,-rpath=$(SRC_DIR)/bindings/Fortran -Wl,-rpath=$(SRC_DIR)/cuda 
LDFLAGS+=-lmiraculix -lfortran 

# CUDA Part
ifeq ($(CUDA), true)
NVCC=$(CUDA_HOME)/bin/nvcc
CXXFLAGS=-Xcompiler --std=c++14 -Xcompiler -fopenmp  -Xcompiler -fPIC -Xcompiler -g -x cu
CXXFLAGS+=-I$(SRC_DIR)/cuda

LDFLAGS+=-L${CUDA_HOME}/lib64 -lcudart -lcublas -lcusolver -lcusparse -lcusolver -ldgemm_compressed_cuda_12
FFLAGS+=-DCUDA
endif

SRC = benchmark.f90
OBJ = $(SRC:.f90=.o)
MODS = $(OBJ:.o=.mod)

LIBFORTRAN = $(SRC_DIR)/bindings/Fortran/libfortran.a
LIBMIRACULIX = $(SRC_DIR)/miraculix/libmiraculix.a

EXE = benchmark.out


all: $(EXE) uint8_lib


$(EXE): $(OBJ) $(LIBMIRACULIX) $(LIBFORTRAN)
	$(FORTRAN) $(FFLAGS) $(OBJ) -o $@ $(LDFLAGS) 

%.o: %.f90
	$(FORTRAN) $(FFLAGS) -c $<

$(OBJ): $(LIBFORTRAN)

uint8_lib: cublas_uint8.o $(SRC_DIR)/cuda/cuda_utils.o
	$(NVCC) -shared -o $(<:.o=.so) $^ -lcublas

%.o: %.cu 
	$(NVCC) $(CXXFLAGS) -c $< -o $@

clean:
	$(RM) $(OBJ) $(MODS) $(EXE) *.optrpt

test: all
	./$(EXE) ../data/small.bed ../data/small.freq
