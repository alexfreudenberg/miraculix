CUDA=true
COMPILER=intel

#Martin's includes:
MK_SRC=./
SRC=./
include makefile.c.mk

ifeq ($(COMPILER), gcc)
 CC=g++
 CCOMPILER=gcc
 CFLAGS=-I/usr/local/include -O2 $(DEBUGGING_FLAG) -mtune=native -march=native
 CCFLAGS=-I/usr/local/include -O2 -fopenmp -pthread $(DEBUGGING_FLAG) -mtune=native -march=native
 LINK_FLAGS=-llapack -lblas -lgfortran -lm -lquadmath -fopenmp

 LINKER=-std=gnu++14 -L/usr/local/lib -L/usr/lib64

else ifeq ($(COMPILER), intel)

 CC=icc
 CFLAGS=-O2 -diag-disable=10441 -fPIC $(DEBUGGING_FLAG) 
 CCFLAGS=-O2 -qopenmp -parallel -pthread -diag-disable=10441 -fPIC $(DEBUGGING_FLAG)
 LINK_FLAGS=-liomp5 -lmkl_blas95_lp64 -lmkl_lapack95_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core  -lpthread -lm -ldl -lstdc++

 LINKER=-std=gnu++14

endif


LIB = libmiraculix.a
SHARED_LIB = miraculix.so

ifeq ($(CUDA), true)
CUDA_LIB = ../cuda/
LINK_FLAGS+=-L${CUDA_HOME}/lib64 -lcudart -lcublas -lcusolver -lcusparse -lcusolver -L$(CUDA_LIB) -ldgemm_compressed_cuda

CCFLAGS+=-DCUDA # Indicates that CUDA functions need to be linked 
FFLAGS+=-DCUDA
LINK=$(CC)
endif



shared_lib: $(SHARED_LIB)
lib: $(LIB)

$(SHARED_LIB): $(ALL)
	$(LINK) -shared $^ -o $@ $(LINK_FLAGS) 

$(LIB): $(ALL)
	ar cr $@ $(ALL)


vtest: all
	valgrind --tool=memcheck --leak-check=full --num-callers=20 ./$(EXE) data/small.bed data/small.freq


vtestt: all
	valgrind --tool=memcheck --leak-check=full --num-callers=20 ./$(EXET) data/small.bed data/small.freq

.PHONY: clean lib
clean: 
	rm *.o