DEBUG=false

CUTLASS_PATH=../cutlass/include/
CUTLASS_UTILS=../cutlass/tools/util/include

cu_nums_11= 70 75 80 
CU_NUMS_11=$(subst $() $(),;,$(cu_nums_11))
CU_FLAGS_11=$(foreach i,$(cu_nums_11),--generate-code arch=compute_$(i),code=sm_$(i)) -DCUTLASS_NVCC_ARCHS='$(CU_NUMS_11)' 

cu_nums_12= 70 75 80 90
CU_NUMS_12=$(subst $() $(),;,$(cu_nums_12))
CU_FLAGS_12=$(foreach i,$(cu_nums_12),--generate-code arch=compute_$(i),code=sm_$(i)) -DCUTLASS_NVCC_ARCHS='$(CU_NUMS_12)' 


NVCC=$(CUDA_HOME)/bin/nvcc
CXXFLAGS=-Xcompiler --std=c++14 -Xcompiler -fopenmp  -Xcompiler -fPIC -Xcompiler -g -x cu -DDO_PARALLEL -I$(CUTLASS_PATH) -I$(CUTLASS_UTILS) -I./ -DCOMMIT_ID=\"$(shell git rev-parse --short HEAD)\" 
LIBS=-L${CUDA_HOME}/lib64 -lcudart -lcublas -lcusolver -lcusparse -lcusolver -Xcompiler -fopenmp

ifeq ($(DEBUG), true)
CXXFLAGS+=-DDEBUG
endif

.PHONY: lib lib12 shared_lib
shared_lib: dgemm_compressed_cuda.o solve_gpu.o cuda_utils.o
	$(NVCC) -shared -o $(<:.o=.so) $^

lib: dgemm_compressed_cuda.o solve_gpu.o cuda_utils.o
	ar rcs lib$(<:.o=.a) $^ 

lib12: dgemm_compressed_cuda_12.o solve_gpu_12.o cuda_utils_12.o
	ar rcs lib$(<:.o=.a) $^ 

%_12.o: %.cu 
	$(NVCC) $(CU_FLAGS_12) $(CXXFLAGS) -c $< -o $@

%.o: %.cu 
	$(NVCC) $(CU_FLAGS_11) $(CXXFLAGS) -c $< -o $@

.PHONY: clean
clean:
	rm *.o *.a