##################################################
## RandomFieldsUtils
##################################################

##CROSS=

MY_SSE2=
MY_SSE3=
MY_SSSE3=
MY_SSE41=
MY_AVX=rfu_avx.o
MY_AVX2=rfu_avx2.o
MY_AVX512F=
MY_CU_FILES=solve_RFU_61.o gpu_info_61.o

MY_CC_FILES=AutoRandomFieldsUtils.o  extern_rfu.o maths.o options_rrfu.o solve_rrfu.o sort.o sortLong.o utils.o win_linux_aux.o xport_import_rrfu.o    

MY_CC_FOR_R_FILES=brdomain_R.o kleinkram.R.o maths.R.o options_RFU.o RFoptions_RFU.o solve_RFU.o sort.R.o utils.R.o xport_import_RFU.o gpu_info_R.o 

MY_CC_DISTRIBUTED=kleinkram.o compatibility.c.o compatibility.sexp.o 
MY_C_FILES=zzz_RFU.o  ## unused

#MY_FORTRAN_FILES=bckslvmodified.f cholmodified.f spamown.f 

