/* 
    Authors
    Alexander Freudenberg, afreuden@mail.uni-mannheim.de

    Copyright (C) 2022-2023  Alexander Freudenberg

    Licensed under the Apache License, Version 2.0 (the "License");
 you may not use this file except in compliance with the License.
 You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

 Unless required by applicable law or agreed to in writing, software
 distributed under the License is distributed on an "AS IS" BASIS,
 WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 See the License for the specific language governing permissions and
 limitations under the License.
*/

#include <unistd.h>
#include <cuda_runtime_api.h>
#include <cuda_runtime.h>
#include <cuda.h>

#include "Basic_utils.h"
#include "errors_messages.h"


#define NUMBER_OF_ATTRIBUTES 5

SEXP gpu_info_61(SEXP DEVICES){
  int N_GPUS = length(DEVICES),
    protects = 0,
    count = 0,
    i = 0,
    device_number = -1,
    *devices = INTEGER(DEVICES);
    
  size_t free_mem = 0,
    total_mem = 0;
  double scaling_constant = 1024.0 * 1024.0 * 1024.0;
  
  cudaDeviceProp deviceProp;
  SEXP gpu_info_list;
  SEXP list_element;
  SEXP dimnames = PROTECT(allocVector(VECSXP, NUMBER_OF_ATTRIBUTES));
  protects++;
  // define dimnames of gpu_info_list
  SET_VECTOR_ELT(dimnames, 0, mkChar("Device Number"));
  SET_VECTOR_ELT(dimnames, 1, mkChar("Device Name"));
  SET_VECTOR_ELT(dimnames, 2, mkChar("Compute Capability"));
  SET_VECTOR_ELT(dimnames, 3, mkChar("Free Memory"));
  SET_VECTOR_ELT(dimnames, 4, mkChar("Total Memory"));

  // check if any device is out of bound
  cudaGetDeviceCount(&count);
  if(count == 0) ERR0("No CUDA devices found.\n");
  for(i=0; i<N_GPUS; i++){
    device_number = devices[i];
    if( device_number < 0 || device_number >= count)
      ERR1("Device out of bound: %d\n", device_number);
  }
  
  // allocate info list
  PROTECT(gpu_info_list = allocVector(VECSXP, N_GPUS));
  protects++;
  
  // Fill list with device info
  for(i = 0; i<N_GPUS; i++){
    PROTECT(list_element = allocVector(VECSXP, NUMBER_OF_ATTRIBUTES));
    
    cudaGetDeviceProperties(&deviceProp, devices[i]);
    cudaSetDevice(devices[i]);
    cudaMemGetInfo(&free_mem, &total_mem);
    
    SET_VECTOR_ELT(list_element, 0, ScalarInteger(devices[i]));
    SET_VECTOR_ELT(list_element, 1, mkString(deviceProp.name));
    SET_VECTOR_ELT(list_element, 2,
		       ScalarInteger(deviceProp.major * 10 +deviceProp.minor));
    SET_VECTOR_ELT(list_element, 3,
		   ScalarReal( ((double)free_mem)/(scaling_constant)));
    SET_VECTOR_ELT(list_element, 4,
		   ScalarReal( ((double)total_mem)/(scaling_constant)));
    
    setAttrib(list_element, R_NamesSymbol, dimnames);
    
    SET_VECTOR_ELT(gpu_info_list, i, list_element);
    UNPROTECT(1);
  }
  
  UNPROTECT(protects);
  return(gpu_info_list);  
}
