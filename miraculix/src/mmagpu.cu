
/*
 Authors 
 Martin Schlather, schlather@math.uni-mannheim.de


 Copyright (C) 2020 -- 2021  Martin Schlather, Alexander Freudenberg

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 3
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.  
*/

#define MY_METHOD MMAGPU
#define BitsPerCode 2 


//////////////////////////////////////////////////
// DO NOT MOVE OR DELETE INCLUDES OR CHANGE ORDER
// very nasty compile errors caused by redefinitions

#include "mmagpuIntern.h"

#include "intrinsics.h"
#include "IntrinsicsBase.h"
#include "error.h"
#include "MX.h"
// #include "options.h"


static void gpuCrossprodIntern(unsigned int* CGM, size_t snps, size_t individuals,
  double* ans, size_t TileSize) {

  // force multiples of 32 byte
  const size_t BytesPerRow = (1 + ((1 + (snps - 1) / CodesPerByte) - 1) / 32) * 32;
  const size_t IntsPerRow = 1 + (BytesPerRow - 1) / sizeof(unsigned int);

  // sanity checks
  // limit Tilesize to individuals
  TileSize = TileSize > individuals ? individuals : TileSize;

  // get total GPU memory
  size_t free_mem;
  cudaMemGetInfo(&free_mem, nullptr);

  // calculates total memory requirements
  size_t req_mem = 2 * BytesPerRow * TileSize * CodesPerByte + TileSize * TileSize * sizeof(unsigned int);
  if(req_mem > free_mem)ERR("Not enough memory available.");

  cudaDeviceProp prop;
  cudaGetDeviceProperties(&prop,0);
  
  PRINTF("Using device %s\n", prop.name);

  // Input data
  data_type* d_x;
  data_type* d_y;
  // Buffer for output
   int32_t* d_val;
  // Buffer for copying back results from device to host
   int32_t* h_val;

   
  const int size_of_input = BytesPerRow * TileSize * CodesPerByte / MEMORY_FACTOR;
  const int size_of_output = sizeof(int32_t) * TileSize * TileSize;
  // Initialization of buffers: We calculate n_streams of tile matrix multiplications in parallel and allocate the corresponding amount of memory
  cudaMalloc((void**)&d_x, n_streams * size_of_input);
  cudaMalloc((void**)&d_y, n_streams * size_of_input);
  cudaMalloc((void**)&d_val, n_streams * size_of_output);
  cudaMallocHost((void**)&h_val, n_streams * size_of_output); 
  err_check("Memory allocation: ");

  //initialization of cutlass gemm operators
  using ColumnMajor = cutlass::layout::ColumnMajor;
  using RowMajor = cutlass::layout::RowMajor;
  using TensorOp = cutlass::arch::OpClassTensorOp;

/* Links:
https://github.com/NVIDIA/cutlass/blob/master/media/docs/functionality.md
https://github.com/NVIDIA/cutlass/blob/master/test/unit/gemm/device/gemm_s8t_s8n_s32n_tensor_op_s32_sm75.cu
    */
    using ElementA_ = data_type;
    using LayoutA_  = RowMajor;
    using ElementB_ = data_type;
    using LayoutB_  = ColumnMajor;
    using ElementC_ = int32_t;
    using LayoutC_  = RowMajor;
    using ElementAccumulator_ = ElementC_;
    using OperatorClass_ = TensorOp;
    using ArchTag_  = cutlass::arch::Sm75;
    using ThreadblockShape_ = typename cutlass::gemm::device::DefaultGemmConfiguration<
    OperatorClass_, ArchTag_, ElementA_, ElementB_, ElementC_,
    ElementAccumulator_>::ThreadblockShape;
    using WarpShape_ = typename cutlass::gemm::device::DefaultGemmConfiguration<
    OperatorClass_, ArchTag_, ElementA_, ElementB_, ElementC_,
    ElementAccumulator_>::WarpShape ;
    using InstructionShape_ = typename cutlass::gemm::device::DefaultGemmConfiguration<
    OperatorClass_, ArchTag_, ElementA_, ElementB_, ElementC_,
    ElementAccumulator_>::InstructionShape;
    using EpilogueOutputOp_ = typename cutlass::gemm::device::DefaultGemmConfiguration<
    OperatorClass_, ArchTag_, ElementA_, ElementB_, ElementC_,
    ElementAccumulator_>::EpilogueOutputOp;
    using ThreadblockSwizzle_ = typename cutlass::gemm::threadblock::GemmIdentityThreadblockSwizzle<>;
    const int Stages = cutlass::gemm::device::DefaultGemmConfiguration<
    OperatorClass_, ArchTag_, ElementA_, ElementB_, ElementC_,
    ElementAccumulator_>::kStages;
    const int AlignmentA = cutlass::gemm::device::DefaultGemmConfiguration<
    OperatorClass_, ArchTag_, ElementA_, ElementB_, ElementC_,
    ElementAccumulator_>::kAlignmentA;
    const int AlignmentB = cutlass::gemm::device::DefaultGemmConfiguration<
    OperatorClass_, ArchTag_, ElementA_, ElementB_, ElementC_,
    ElementAccumulator_>::kAlignmentB;
    const bool SplitKSerial = false;
    using Operator_ = typename cutlass::gemm::device::DefaultGemmConfiguration< OperatorClass_, ArchTag_, ElementA_, ElementB_, ElementC_,    ElementAccumulator_>::Operator;
    const bool IsBetaZero = false;

    /* Links:
    https://github.com/NVIDIA/cutlass/blob/master/media/docs/functionality.md
    https://github.com/NVIDIA/cutlass/blob/master/test/unit/gemm/device/gemm_s8t_s8n_s32n_tensor_op_s32_sm75.cu
    */


    using CutlassGemm = cutlass::gemm::device::Gemm<
        ElementA_,        // Data-type of A matrix
        LayoutA_,         // Layout of A matrix
        ElementB_,        // Data-type of B matrix
        LayoutB_,      // Layout of B matrix
        ElementC_,              // Data-type of C matrix
        LayoutC_,        // Layout of C matrix
        ElementAccumulator_,
        OperatorClass_,            // tag indicating Tensor Cores,          // tag indicating Tensor Cores
        ArchTag_,
        ThreadblockShape_,
        WarpShape_,
        InstructionShape_,
        EpilogueOutputOp_ ,
        ThreadblockSwizzle_,
        Stages,
        AlignmentA,
        AlignmentB,
        SplitKSerial,
        cutlass::arch::CustomOp,
        IsBetaZero  >;

  // Define a CUTLASS GEMM type
  CutlassGemm gemm_operator;
  
  // Main loop
  // Calculates matrix multiplications in parallel: Each thread in this loop sends its data to a different stream on the device. The threads calculate concurrently and send the output back to main memory. Memory copies are asynchronous to take full advantage of the memory bandwidth.
  #ifdef DO_PARALLEL
  #pragma omp parallel for num_threads(n_streams <  1+(individuals-1)/TileSize ? n_streams : 1+(individuals-1)/TileSize) schedule(dynamic)   
  #endif
  for (int64_t i = 0; i < individuals; i += TileSize) {
      int threadidx = omp_get_thread_num();
      cudaStream_t stream;
      cudaError_t err = cudaStreamCreate ( & stream );
      cudaStreamSynchronize(stream);

      if(err != cudaSuccess)
          ERR("Stream couldn't be created");
      
      // Pointer to the first element of current rows
      unsigned int *x = (CGM + i * IntsPerRow);
      data_type *x_dev = d_x + threadidx * TileSize * BytesPerRow ;
      data_type *y_dev = d_y + threadidx * TileSize * BytesPerRow ;

      // Number of rows in matrix
      size_t const rows_left = individuals - i;
      // Size x of current tile
      size_t const x_tile_size = TileSize < rows_left ? TileSize : rows_left;

    // Legacy code for 8-bit multiplication
     /* thrust::transform(thrust::cuda::par.on(stream), reinterpret_cast<unsigned char*>(x), reinterpret_cast<unsigned char*>(x) + x_tile_size * BytesPerRow,
          thrust::cuda::pointer<dev_data_type>(reinterpret_cast<dev_data_type*>(raw_x + threadidx * TileSize * BytesPerRow * CodesPerByte / memory_constant )), unpack_data<data_type>());
      */    
     cudaMemcpyAsync(x_dev, x, x_tile_size * BytesPerRow,cudaMemcpyHostToDevice, stream);

          cudaStreamSynchronize(stream);
          err_check("Copy 1:");

      // Inner loop
      for (int64_t j = i; j < individuals; j += TileSize) {

          // Same as above with y
          size_t const columns_left = individuals - j;
          size_t const y_tile_size = TileSize < columns_left ? TileSize : columns_left;
          unsigned int* y = (CGM + j * IntsPerRow);

    // Legacy code for 8-bit multiplication
     /*     thrust::transform(thrust::cuda::par.on(stream), reinterpret_cast<unsigned char*>(y), reinterpret_cast<unsigned char*>(y) + y_tile_size * BytesPerRow,
              thrust::cuda::pointer<dev_data_type>(reinterpret_cast<dev_data_type*>(raw_y + threadidx *  TileSize * BytesPerRow * CodesPerByte /memory_constant)), unpack_data<data_type>());
      */   
      
         cudaMemcpyAsync(y_dev, y, y_tile_size * BytesPerRow, cudaMemcpyHostToDevice, stream);
          err_check("Copy 2:");
          cudaStreamSynchronize(stream);

          //initialize gemm arguments
          CutlassGemm::Arguments args(
              { int(x_tile_size), int(y_tile_size), int(BytesPerRow * CodesPerByte  / COMPRESSION_GPU) },
              { x_dev, int(BytesPerRow * CodesPerByte / COMPRESSION_GPU) }, 
              { y_dev, int(BytesPerRow * CodesPerByte  / COMPRESSION_GPU) },
              { d_val + threadidx * TileSize * TileSize, int(y_tile_size) },  
              { d_val + threadidx * TileSize * TileSize, int(y_tile_size) },  
              { 1, 0 });                     
          cudaStreamSynchronize(stream);

          // compute Multiplication
          cutlass::Status status = gemm_operator(args, nullptr, stream);
          cudaStreamSynchronize(stream);
          err_check("Calculation:");


          // Copy results back to host
          cudaMemcpyAsync(h_val + threadidx * TileSize * TileSize, d_val + threadidx * TileSize * TileSize, TileSize * TileSize * sizeof(int32_t), cudaMemcpyDeviceToHost, stream);
          err_check("Copying back:");
          
          cudaStreamSynchronize(stream);

          if(*(h_val + threadidx * TileSize * TileSize)==0){printf("Computation failed at thread %d, (%d,%d)\n",threadidx, i,j);
                      print_kernel <<<1,1>>> ((int32_t*) d_val + threadidx * TileSize * TileSize);
              j-=TileSize;
              continue;
          }
          err_check("Copy back:");

          // Loop over tile and store values in output matrix 
          #ifdef DO_PARALLEL
          #pragma omp parallel for num_threads(n_streams) schedule(static)   
          #endif
          for (int64_t di = 0; di < x_tile_size; ++di)
          {
              for (int64_t dj = 0; dj < y_tile_size; ++dj)
              {
                  // Get result
                  const auto Mij = *(h_val +threadidx * TileSize * TileSize + dj + di * y_tile_size ); 
                  // Create pointers to the output matrix (because it is symmetric we use two pointers)
                  double* ans0 = ans + (i + di),
                      * ans1 = ans + (i + di) * individuals;
                  // Store result in ouput matrix
                  ans0[(j + dj) * individuals] = ans1[j + dj] = (double) Mij;
              }
             // PRINTF("\n");
          }
 
      }
      cudaStreamDestroy(stream);
  }

  // Free memory
  cudaFree(d_x);
  cudaFree(d_y);
  cudaFree(d_val);
  cudaFreeHost(h_val);

}


static SEXP matrix_start_Intern(Uint snps, Uint individuals, SEXP VARIABLE_IS_NOT_USED file){
SEXP Code;
PROTECT(Code = CreateEmptyCodeVector(snps, individuals, MY_METHOD));
UNPROTECT(1);
return Code;
}


static void crossprodIntern(Uint *CM, Uint snps, Uint individuals,
        double *ans) {
// tilse_size needs to stay the same: for smaller values we experience undocumented calculation failures on the device
const size_t tilesize = 2048; 

// Initialize host pointers and copy input data cuda managed memory
Uint* h_CM; 
const size_t BytesPerIndiv = UnitsPerIndiv256(snps) * BytesPerUnit;  cudaMallocHost((void**)&h_CM, individuals * BytesPerIndiv);
MEMCOPY(h_CM, CM, individuals * BytesPerIndiv);

gpuCrossprodIntern(h_CM, snps, individuals, ans, tilesize);
cudaFreeHost(h_CM);

}


  // Test for right compute capability
  bool check_7_5() {
    cudaDeviceProp deviceProp;  
    cudaGetDeviceProperties(&deviceProp, 0);
    if((deviceProp.major < 7) | (deviceProp.major == 7 & deviceProp.minor <5))
    {  
      ERR("No GPU kernel for compute capability 7.0 or lower yet.");
      return false;
    }
    return true;
        // helpful values
  }
  

void crossprod_mmagpu(Uint *CGM, Uint snps, Uint individuals, double *ans){
  crossprodIntern(CGM, snps, individuals, ans);
}

SEXP matrix_start_mmagpu( Uint snps,Uint individuals, SEXP file) {
  check_7_5();
  return matrix_start_Intern(snps, individuals, file);
}
 
bool useMMAGPU(snpcoding method) {
  return method == MMAGPU ? true : false;  
}




void gpu_relmat_cublas(Uint* M, double* A, Uint snps, Uint individuals){
    /*
        Calculates the crossproduct of M with cublas and stores the result in A.
        Input:
            M: non-encoded matrix of dimension snps x indiv (k x n) storing genomic information
            A: pointer to result matrix
            snps: Number of snps
            individuals: number of individuals
        Output:
            A: matrix containing the type-casted result of M^T * M
        
        Note: cublas is fortran based and therefore assumes M is column-major. Therefore to calculate
            A we instruct cublasgemmex to calculate M * M^T and adjust its parameters.
            Furthermore, gemmex requires the matrix M to have a row number that is a multiple of four
            Therefore this function implements a zero-padding to add extra rows
    */
    
    //Define auxiliary variables as needed for gemmex
        Uint n = individuals;
        Uint m = individuals;
        Uint k = snps;
    
    //Auxiliary padding variables for padding
        Uint k_pad_diff = (PADDIM - k % PADDIM) % PADDIM;
        Uint k_pad = k + k_pad_diff;
        Uint dim = m * k_pad;
    
    //Start timing copy and calculation time
    #ifdef DEBUG
        std::chrono::time_point<std::chrono::high_resolution_clock> timer_start;
        std::chrono::time_point<std::chrono::high_resolution_clock> timer_stop;
        timer_start = std::chrono::high_resolution_clock::now();
    #endif
    //Declare cublas variables and allocate memory
        cublasHandle_t handle;
        cublasCreate(&handle);
        int8_t *d_M, *h_M;
        int32_t *d_C, *h_C;
        int32_t alpha = 1.f;
        int32_t beta = 0.f;
        cudaMalloc(&d_M, sizeof(int8_t) * dim);
        cudaMalloc(&d_C, sizeof(int32_t) * n * m );
        cudaMallocHost((void **)&h_M, sizeof(int8_t) * dim);
        cudaMallocHost((void **)&h_C, sizeof(int32_t) * n * m);
    
    
    
    //Type-cast matrix M to int8 and store the result in page-locked memory
    //Zero-pad matrix to get a row number that is a multiple of four
    #ifdef DO_PARALLEL
    #pragma omp parallel for num_threads(CORES)   
    #endif
        for(int i = 0; i < n; i++){
            for(int j = 0; j < k_pad; j++){
            h_M[j + i * k_pad] = (int8_t) (j< k ?  M[j + i * k] : 0 );
            }
        }
    
    
    //Copy int8 matrix to device
    cudaMemcpy(d_M, h_M, sizeof(int8_t) * dim, cudaMemcpyHostToDevice);  

    //Calculate the crossproduct and check for errros
        cublasStatus_t stat = cublasGemmEx(handle,
            CUBLAS_OP_T,
            CUBLAS_OP_N,
            n,
            m,
            k_pad,
            &alpha,
            d_M,
            CUDA_R_8I,
            k_pad, // I have no idea why this doesnt need to be individuals, same below
            d_M,
            CUDA_R_8I,
            k_pad,
            &beta,
            d_C,
            CUDA_R_32I,
            n,
            CUDA_R_32I, //CUBLAS_COMPUTE_32I,
            CUBLAS_GEMM_DEFAULT
            );
        
        if(stat) PRINTF("GemmEx failed.");
        cudaDeviceSynchronize();
    
    
    //copy result back to host
        cudaMemcpy(h_C, d_C, sizeof(int32_t) * n * m, cudaMemcpyDeviceToHost);
    
    //Convert result to double and store it in output matrix A
    #ifdef DO_PARALLEL
    #pragma omp parallel for num_threads(CORES)   
    #endif
        for (int i = 0; i < n * m; i++) A[i] = (double) h_C[i];
    
    //Free memory 
        cublasDestroy(handle);
        cudaFree(d_M);
        cudaFree(d_C);
        cudaFreeHost(h_C);
        cudaFreeHost(h_M);
    
    //Stop timer
    #ifdef DEBUG
        timer_stop = std::chrono::high_resolution_clock::now();
        PRINTF("Time: %.3f s\n", ((float) std::chrono::duration_cast<std::chrono::microseconds>(timer_stop - timer_start).count())/1000000.0 );
    #endif
    } 
/*
This is a possible alternative implementation in cublasLt. Might be faster?
    cublasLtCreate(&handle);

    cublasLtMatmulDescCreate(&operationDesc, CUBLAS_COMPUTE_32I, CUDA_R_32I);
    cublasLtMatmulDescSetAttribute(operationDesc, CUBLASLT_MATMUL_DESC_TRANSB, &transb, sizeof(transb));
    cublasLtMatmulDescSetAttribute(operationDesc, CUBLASLT_MATMUL_DESC_TRANSB, &transb, sizeof(transb));

    cublasLtMatrixLayoutCreate(&Adesc, CUDA_R_8I, snps, individuals, snps);
    cublasLtMatrixLayoutCreate(&Bdesc, CUDA_R_8I, individuals, snps, snps);
    cublasLtMatrixLayoutCreate(&Cdesc, CUDA_R_32I, individuals, individuals, individuals);


    cublasLtMatmulPreferenceCreate(&preference);
    cublasLtMatmulPreferenceInit(preference);

//Hopefully not needed:
    cublasLtMatmulPreferenceSetAttribute(preference, CUBLASLT_MATMUL_PREF_MAX_WORKSPACE_BYTES, &workspaceSize, sizeof(workspaceSize));

    cublasLtMatmulAlgoGetHeuristic(handle, operationDesc, Adesc, Bdesc, Cdesc, Cdesc, preference, 1, &heuristicResult, &returnedResults);

    if (returnedResults == 0) {
        ERR1("Status %s", CUBLAS_STATUS_NOT_SUPPORTED);
    }

    cublasLtMatmul(handle,
        operationDesc,
        d_alpha,  // 1
        d_M,
        Adesc,
        d_M,
        Bdesc,
        d_beta,   // 0
        d_C,
        Cdesc,
        d_C,
        Cdesc,
        &heuristicResult.algo,
        workspace,
        workspaceSize,
        0);

    cublasLtMatmulDescDestroy(operationDesc);
    cublasLtDestroy(handle);
*/