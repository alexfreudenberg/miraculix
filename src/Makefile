CUDA=true

TARGETS=libmiraculix shared_lib libfortran

ifeq ($(CUDA), true)
TARGETS+=libmiraculixGPU
endif


.PHONY: libmiraculixGPU libmiraculix shared_lib libfortran check_cuda clean
all: $(TARGETS)

shared_lib: libmiraculixGPU
	cp Makefile_miraculix ./miraculix/Makefile
	cd miraculix && make shared_lib

libmiraculix:
	cp Makefile_miraculix ./miraculix/Makefile
	cd miraculix && make lib

libmiraculixGPU: check_cuda
	cd cuda && make lib

libfortran: libmiraculix libmiraculixGPU
	cd bindings/Fortran && make lib

check_cuda: 
ifndef CUDA_HOME
$(error Environment variable CUDA_HOME not set - check your CUDA installation.)
endif
CUDA_FILES=bin/nvcc lib64/libcudart.so lib64/libcublas.so lib64/libcusolver.so lib64/libcusparse.so
$(foreach file,$(CUDA_FILES),$(if $(wildcard $(CUDA_HOME)/$(file)),,$(error CUDA file $(CUDA_HOME)/$(file) does not exist)))


clean:
	$(RM) cuda/*.o cuda/*.a
	$(RM) miraculix/*.o miraculix/*.a miraculix/*.so
