
# compilation flags
ifeq ($(COMPILER), gcc)
CC=g++
Ccompiler=gcc
FORTRAN=gfortran
CFLAGS=-I/usr/local/include -O2 $(DEBUGGING_FLAG) 
CCFLAGS=-I/usr/local/include -O2 -fopenmp -pthread $(DEBUGGING_FLAG)
FFLAGS=-O2 
LINK_FLAGS=-llapack -lblas -lgfortran -lm -lquadmath -fopenmp




SSSE3_FLAG=-mssse3 
SSE41_FLAG=-msse4.1 
AVX_FLAG=-mavx
#AVX2_FLAG=-mavx2 
AVX2_FLAG=-mavx2
AVX512_FLAG=-mavx512f
MAX_AVX512_FLAG=-mavx512f
LINKER=-std=gnu++14 -L/usr/local/lib -L/usr/lib64

else ifeq ($(COMPILER), intel)

CC=icc
Ccompiler=icc
FORTRAN=ifort
CFLAGS=-O2 -diag-disable=10441 $(DEBUGGING_FLAG) 
CCFLAGS=-O2 -qopenmp -parallel -pthread -diag-disable=10441 $(DEBUGGING_FLAG)
FFLAGS=-O2 
FFLAGS+=-I${MKLROOT}/include -I${MKLROOT}/include/intel64/lp64 -L${MKLROOT}/lib/intel64
LINK_FLAGS=-liomp5 -lmkl_blas95_lp64 -lmkl_lapack95_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core  -lpthread -lm -ldl

SSSE3_FLAG=-xSSSE3
SSE41_FLAG=-xSSE4.1
AVX_FLAG=-xAVX
AVX2_FLAG=-xCORE-AVX2
AVX512_FLAG=-xCOMMON-AVX512
MAX_AVX512_FLAG=-xCORE-AVX512
LINKER=-std=gnu++14

else ifeq ($(COMPILER), clang)

CC=clang++
Ccompiler=clang
FORTRAN=gfortran
SANITIZE_FLAGS=-fsanitize=undefined -fsanitize=integer -fsanitize=address -fno-sanitize=float-divide-by-zero -fno-sanitize=alignment -fno-omit-frame-pointer -frtti
SANITIZE_FLAGS=-fsanitize=address -fno-sanitize=alignment
CFLAGS=-I/usr/local/include -O2 $(DEBUGGING_FLAG) $(SANITIZE_FLAGS)
CCFLAGS=-I/usr/local/include -O2 -fopenmp -pthread $(DEBUGGING_FLAG) $(SANITIZE_FLAGS)
FFLAGS=-O2 
LINK_FLAGS=-llapack -lblas -lgfortran -lm -lquadmath -fopenmp -lasan -L /usr/local/lib 

CCFLAGS+   format -Wformat-security -Wformat-y2k -Wuninitialized -Wall -Wextra -Wshadow -Wpointer-arith -pedantic -Wswitch-default -march=native -Wno-unused-variable -Wno-unused-function -mtune=native -Wno-ignored-attributes -Wno-deprecated-declarations -Wno-parentheses  -Wformat-overflow=2

CCFLAGS+=-g  -pipe


SSSE3_FLAG=-mssse3  
AVX_FLAG=-mavx
#AVX2_FLAG=-mavx2 
AVX2_FLAG=-mavx2
AVX512_FLAG=-mavx512 
MAX_AVX512_FLAG=-mavx512
LINKER=-std=gnu++14 -L/usr/local/lib -L/usr/lib64


endif




# Flags for internal purposes
PURPOSE_FLAG=-DREQUIRED_SIMD=2 -DLINUXCPUID

include $(MK_SRC)/makefile.rfu.mk
RFU_SSSE3:=$(MY_SSSE3)
RFU_SSE41:=$(MY_SSE41)
RFU_AVX:=$(MY_AVX) 
RFU_AVX2:=$(MY_AVX2)
RFU_AVX512F:=$(MY_AVX512F)
RFU_MAX_AVX512:=$(MY_MAX_AVX512F)
RFU_CU:=$(MY_CU)
RFU_CC:=$(MY_CC_FILES)
RFU_CC_R:=$(MY_CC_FOR_R_FILES)
RFU_F:=$(MY_FORTAN_FILES)


include $(MK_SRC)/makefile.miraculix.mk
ALL_SSSE3:=$(MY_SSSE3) $(RFU_SSSE3)
ALL_SSE41:=$(MY_SSE41) $(RFU_SSE41)
ALL_AVX:=$(MY_AVX) $(RFU_AVX) 
ALL_AVX2:=$(MY_AVX2) $(RFU_AVX2)
ALL_AVX512F:=$(MY_AVX512F) $(RFU_AVX512F)
ALL_MAX_AVX512F:=$(MY_MAX_AVX512F) $(RFU_MAX_AVX512F)
ALL_CU:=$(MY_CU) $(RFU_CU)
ALL_CC:=$(MY_CC_FILES) $(RFU_CC)
ALL_CC_R:=$(MY_CC_FOR_R_FILES) $(RFU_CC_R)
ALL_F:=$(MY_FPRTRAN_FILES) $(RFU_F)
ALL_C:=$(MY_C_FILES)

ALL:=$(ALL_CC) $(ALL_SSSE3) $(ALL_SSE41) $(ALL_AVX) $(ALL_AVX2)  $(ALL_AVX512F) $(ALL_MAX_AVX512F) $(ALL_F) $(ALL_C)


alle : msg $(ALL)

allparallel : msg_parallel msg do_parallel Main

msg_parallel :
	echo ******** PARALLEL ******** 

msg :
	@echo ""
	@echo SRC=$(SRC)


intrinsics.o: $(SRC)/intrinsics.cc
	$(CC) -c $< -o $@ $(CCFLAGS) $(AVX2_FLAG)

%avx2.o: $(SRC)/%avx2.cc
	$(CC) -c $< -o $@ $(CCFLAGS) $(AVX2_FLAG)

%avx.o: $(SRC)/%avx.cc
	$(CC) -c $< -o $@ $(CCFLAGS) $(AVX_FLAG)

#%512.o: $(SRC)/%512.cc
#	$(CC) -c $< -o $@ $(CCFLAGS) $(AVX512_FLAG) 

%256.o: $(SRC)/%256.cc
	$(CC) -c $< -o $@ $(CCFLAGS) $(AVX2_FLAG)

%128.o: $(SRC)/%128.cc
	$(CC) -c $< -o $@ $(CCFLAGS) $(SSSE3_FLAG)

%128_41.o: $(SRC)/%128_41.cc
	$(CC) -c $< -o $@ $(CCFLAGS) $(SSE41_FLAG)

%.o: $(SRC)/%.cc
	$(CC) -c $< -o $@ $(CCFLAGS)

%.o: $(SRC)/%.c
	$(Ccompiler) -c $< -o $@ $(CFLAGS)

%.o: $(NL)/%.f90
	$(FORTRAN) $(FFLAGS) -c $< -o $@


%avx.o: $(RFU)/%avx.cc
	 $(CC) -c $< -o $@ $(CCFLAGS) $(AVX_FLAG)

%avx2.o: $(RFU)/%avx2.cc
	 $(CC) -c $< -o $@ $(CCFLAGS) $(AVX2_FLAG)

%.o: $(RFU)/%.cc
	$(CC) -c $< -o $@ $(CCFLAGS)

%.o: $(RFU)/%.c
	$(Ccompiler) -c $< -o $@ $(CFLAGS)

%.o: $(RFU)/%.f
	$(FORTRAN) -c $< -o $@ $(FFLAGS)


Main: 
	$(CC) $(LINKER) -o run_$(COMPILER) $(ALL) $(LINK_FLAGS)


cleaner:
	rm *.o


ALL_CC_R_CC=$(ALL_CC_R:%.o=%.cc)
rmR : 
	echo $(ALL_CC_R_CC)
	rm  $(ALL_CC_R_CC) -f

ALL_CC_R_CC:
	echo $(ALL_CC_R)
	echo $*.*
#	rm $*.*
