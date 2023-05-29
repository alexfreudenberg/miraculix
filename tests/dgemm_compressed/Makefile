
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

SRC_DIR=$(PWD)/../../src
DATA_DIR=../../data

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

EXE = test_5codesapi.out
EXET = test_5codesapi_t.out



all: $(EXE) $(EXET)

$(LIBMIRACULIX):
	@cd $(SRC_DIR) && make libmiraculix

$(EXE): $(EXE:.out=.o) $(LIBMIRACULIX) $(LIBFORTRAN)
	$(FORTRAN) $(FFLAGS) $< -o $@ $(LDFLAGS) 

$(EXET): $(EXE:.out=.o) $(LIBMIRACULIX) $(LIBFORTRAN)
	$(FORTRAN) $(FFLAGS) $< -o $@ $(LDFLAGS)

%.o: %.f90
	$(FORTRAN) $(FFLAGS) -c $<

$(EXE).o: $(LIBFORTRAN)
$(EXET).o: $(LIBFORTRAN)

clean:
	$(RM) $(OBJ) $(MODS) $(LIB) $(EXE) $(EXET) $(ALL) *.optrpt

tests: all
	./$(EXE) $(DATA_DIR)/small.bed $(DATA_DIR)/small.freq

testm: all
	./$(EXE) $(DATA_DIR)/medium.bed $(DATA_DIR)/medium.freq


testl: all
	./$(EXE) $(DATA_DIR)/large.bed $(DATA_DIR)/large.freq


testt: all  
	./$(EXET) $(DATA_DIR)/small.bed $(DATA_DIR)/small.freq


vtest: all
	valgrind --tool=memcheck --leak-check=full --num-callers=20 ./$(EXE) data/small.bed data/small.freq


vtestt: all
	valgrind --tool=memcheck --leak-check=full --num-callers=20 ./$(EXET) data/small.bed data/small.freq
