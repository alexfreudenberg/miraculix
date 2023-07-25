CUDA=false
COMPILER=intel

#miraculix source files:
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

 CC=icx
 CFLAGS=-O2 -fPIC $(DEBUGGING_FLAG) -fiopenmp
 CCFLAGS=-O2 -qopenmp -pthread -fPIC $(DEBUGGING_FLAG) -fiopenmp
 LINK_FLAGS=-liomp5 -lmkl_blas95_lp64 -lmkl_lapack95_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core  -lpthread -lm -ldl -lstdc++

 LINKER=-std=gnu++14

endif


LIB = libmiraculix.a
SHARED_LIB = miraculix.so

ifeq ($(CUDA), true)
CUDA_LIB = ../cuda/
LINK_FLAGS+=-L${CUDA_HOME}/lib64 -lcudart -lcublas -lcusolver -lcusparse -lcusolver 
# To expose CUDA functions which aren't use by core miraculix, we need to link the library as a whole archive
LINK_FLAGS+=-L$(CUDA_LIB) -Wl,--whole-archive -ldgemm_compressed_cuda -Wl,--no-whole-archive 


CCFLAGS+=-DCUDA
FFLAGS+=-DCUDA
LINK=$(CC)
endif


shared_lib: $(SHARED_LIB)
lib: $(LIB)

$(SHARED_LIB): $(ALL)
	$(LINK) -shared $^ -o $@ $(LINK_FLAGS) 

$(LIB): $(ALL)
	ar cr $@ $(ALL)

clean: 
	rm *.o *.a *.so

.PHONY: clean lib shared_lib
