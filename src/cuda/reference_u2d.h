#include <algorithm>
#include <cctype>
#include <iomanip>
#include <iostream>
#include <random>
#include <string>
#include <vector>

#ifdef MKL
#include "mkl.h"
#include "mkl_blas.h"
#include "omp.h"
#endif

// definition of a type which packs 4 doubles
struct PackedDouble {
    double entries[4];
};

// template <
//     typename ElementA,
//     typename LayoutA,
//     typename ElementB,
//     typename LayoutB,
//     typename ElementC,
//     typename LayoutC>
// void Reference_Gemm(TensorView<ElementA, LayoutA> host_A, TensorView<ElementB, LayoutB> host_B, TensorView<ElementC, LayoutC> host_C);

template <typename Element>
struct Reference_Gemm;

// Reference implementation of u2 x F64 gemm
// u2 data is casted to double beforehand
template <>
struct Reference_Gemm<double> {
    long n, m, k;
    Reference_Gemm(
        long m_,
        long n_,
        long k_)
        : m(m_)
        , n(n_)
        , k(k_)
    {
    }

    // Operator uses unpacked data in A and B
    void operator()(
        double* result,
        uint8_t* host_A_full,
        double* host_B_full)
    {
#ifdef MKL
        printf("Allocating %.1lf GB\n", m * k * sizeof(double)/pow(1024,3.0));
        std::vector<double> host_A_double(m * k);

        printf("Typecast\n");
#pragma omp parallel num_threads(omp_get_max_threads())
        for (long i = 0; i < m; i++) {
          for (int j = 0; j < k; j++) {
            host_A_double[i + j * m] = double(host_A_full[i * k + j]);
          }
        }
        printf("Calculate\n");
        cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, m, n, k, 1.0, host_A_double.data(), m, host_B_full, k, 0, result, m);
#else
        assert(0);
#endif
    }
};

// Reference implementation of u2 x F64 gemm
// Data format: 4 * u2 packed in u8, 4 * F64 packed in PackedDouble
template <>
struct Reference_Gemm<PackedDouble> {
    int n, m, k;
    Reference_Gemm(
        long n_,
        long m_,
        long k_)
        : n(n_)
        , m(m_)
        , k(k_)
    {
    }

    void operator()(
        double* result,
        uint8_t* host_A,
        PackedDouble* host_B)
    {
        for (long i = 0; i < m; i++) {
            for (long j = 0; j < n; j++) {
                for (long l = 0; l < k; l++) {
                    double tmp_res = 0.0;
                    uint8_t tmp_a = host_A[i * k + l];
                    PackedDouble tmp_b = host_B[l + j * k];
                    // printf("%.1lf %.1lf %.1lf %.1lf\n", tmp_b.entries[0],tmp_b.entries[1], tmp_b.entries[2], tmp_b.entries[3]);

                    uint8_t tmp = tmp_a & 0x03;
                    tmp_res += ((tmp == 0x02) ? 2.0 : ((tmp == 0x01) ? 1.0 : 0.0)) * tmp_b.entries[0];

                    tmp = ((tmp_a >> 2) & 0x03);
                    tmp_res += ((tmp == 0x02) ? 2.0 : ((tmp == 0x01) ? 1.0 : 0.0)) * tmp_b.entries[1];

                    tmp = ((tmp_a >> 4) & 0x03);
                    tmp_res += ((tmp == 0x02) ? 2.0 : ((tmp == 0x01) ? 1.0 : 0.0)) * tmp_b.entries[2];

                    tmp = ((tmp_a >> 6) & 0x03);
                    tmp_res += ((tmp == 0x02) ? 2.0 : ((tmp == 0x01) ? 1.0 : 0.0)) * tmp_b.entries[3];

                    result[j * n + j] += tmp_res;
                }
            }
        }
    }
};
