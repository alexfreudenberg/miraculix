#include <algorithm>
#include <cctype>
#include <iomanip>
#include <iostream>
#include <random>
#include <string>
#include <vector>

#include "reference_u2d.h"

int main()
{
    typedef uint8_t d_type;
    int n = 2048; // length of vector
    int length_a = n * sizeof(d_type) * 8 / 2; // length of packed {0,1,2} vector

    // Initialize random devices to simulate {0,1,2} and double vector
    std::random_device rd {};
    std::mt19937 gen { rd() };
    std::normal_distribution<> rnorm { 0, 2 };
    std::uniform_int_distribution<unsigned short> runif(0, 2);

    // Initialize vectors to be filled
    double result1 = 0.0,
           result2 = 0.0;
    std::vector<double> result_vector(n / 4); // Vector to be filled with intermediate results
    std::vector<double> data_double(n); // double vector
    std::vector<uint8_t> a_full(n); // Initial unpacked vector of {0,1,2} to validate scalar product

    // Simulate random numbers for both vectors
    for (int i = 0; i < n; i++)
        data_double[i] = rnorm(gen);
    for (int i = 0; i < n; i++)
        a_full[i] = (uint8_t)runif(gen);

    // Slow but safe version of scalar product
    result1 = std::inner_product(data_double.begin(), data_double.end(), a_full.begin(), 0.0);
    printf("Result1: %lf\n", result1);

    // Convert arrays to std::vectors
    PackedDouble* ptr_pdouble = reinterpret_cast<PackedDouble*>(data_double.data());
    std::vector<PackedDouble> data(ptr_pdouble, ptr_pdouble + n / 4);
    std::vector<uint8_t> a8(n / 4);

    // Pack { 0, 1, 2 } vector into 2 - bits stored in a8
    std::fill(a8.begin(), a8.end(), 0x00);
    for (int i = 0; i < n; i++) {
        a8[i / 4] |= (a_full[i] & 0x03) << ((i % 4) * 2);
    }
    // Calculate scalar product without unpacking a8
    std::transform(data.begin(), data.end(), a8.begin(), result_vector.begin(), [](PackedDouble d, uint8_t a) {
        double res = 0.0;
        // Unpack elements of a8 recursively and convert its elements into doubles
        uint8_t factor = (a & 0x03);
        res += ((factor == 0x02) ? 2.0 : ((factor == 0x01) ? 1.0 : 0.0)) * d.entries[0];
        factor = ((a >> 2) & 0x03);
        res += ((factor == 0x02) ? 2.0 : ((factor == 0x01) ? 1.0 : 0.0)) * d.entries[1];
        factor = ((a >> 4) & 0x03);
        res += ((factor == 0x02) ? 2.0 : ((factor == 0x01) ? 1.0 : 0.0)) * d.entries[2];
        factor = ((a >> 6) & 0x03);
        res += ((factor == 0x02) ? 2.0 : ((factor == 0x01) ? 1.0 : 0.0)) * d.entries[3];

        return res;
    });

    // Accumulate scalar product into result2
    result2 = std::accumulate(result_vector.begin(), result_vector.end(), 0.0, std::plus<double>());
    printf("Result2 %lf\n", result2);

    double result3 = 0.0;
    Reference_Gemm ref_gemm(1, 1, n / 4);
    ref_gemm(&result3, a8.data(), data.data());
    printf("Result3 %lf\n", result3);
    
}