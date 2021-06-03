# r_cuda_skeleton
This is an updated version of rpud by Chi Yau, see https://github.com/cran/rpud.

**Major changes include:**
* Updated registration: C and CUDA functions are now registered as native routines with R. See ["Writing R Extensions"](https://cran.r-project.org/doc/manuals/R-exts.html#Registering-native-routines) for advantages.
* Added support for compilation of C functions outside of CUDA-files: Compilation of .c and .cu is done separately, object files are linked by nvcc.
* Fixed linking bugs: Newer versions of R require more linking flags which couldn't be processed by the original Makefile

**To Dos**
* Try removing GNU flavors from Makefiles to comply with CRAN policies
* Add support for different compilers
* Enable processing of compilation/linking flags 
* Write instructions, examples, documentation, etc.
 
