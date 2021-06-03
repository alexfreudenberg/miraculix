#!/bin/bash

## This script can be used for profiling the CUDA part of miraculix. Replace TEST.R by relevant file
R CMD INSTALL RandomFieldsUtils --configure-args="USE_GPU=yes"
R CMD INSTALL miraculix --configure-args="CXX_FLAGS='-mavx2 -DGPU_DEV' USE_GPU='yes'"
nvprof --profile-child-processes -o test%p.prof 
Rscript TEST.R  

