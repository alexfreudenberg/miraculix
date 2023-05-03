##################################################
## miraculix
##################################################

MY_SSE2=
MY_SSE3=
MY_SSSE3=1bit128.o 2bit128.o
MY_SSE41=5codes128_41.o
MY_AVX=
MY_AVX2=5codes256.o  Scalar256.o 1bit256.o  2bit256.o plink256.o
MY_AVX512F=5codes512.o 1bit512.o 2bit512.o 
MY_CU_FILES=mma_75.o mma_61.o solve_61.o
MY_MAX_AVX512F=intrinsics.o

MY_CC_FILES=plinkUint.o 5codesChar.o 5codesUint.o kleinkram.o FilesUint.o Files32.o transformUint.o Haplo2_Uint.o HaploUint.o OneByteUint.o OneByte32.o OneByte64.o extern.o  haplogeno.o  mmagpuDummy.o 4ByteUint.o Automiraculix.o options.o  1bitUint.o 1bit32.o  1bit64.o 2bitUint.o 2bit32.o 2bit64.o  3bitUint.o 3bit64.o Scalar.o MX.o  xport_import.o compatibility.c.o compatibility.sexp.o  Vector.matrix.Uint.o Vector.matrix.LD.o Vector.matrix.D.o 

MY_CC_FOR_R_FILES=transformR.o kleinkram.R.o HaploR.o haplogeno.R.o Vector.matrix.R.o options.R.o ScalarR.o utils_miraculix.R.o windower.R.o MoBPS_R.o xport_import.R.o scan.R.o

MY_CC_DISTRIBUTED=
MY_C_FILES=5codesAPI.o zzzR.o 

MY_FORTRAN_FILES=