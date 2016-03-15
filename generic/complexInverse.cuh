#ifndef COMPLEX_INVERSE_CUH
#define COMPLEX_INVERSE_CUH

#include <cuComplex.h>

__device__ void getComplexLU (const int, const int, cuDoubleComplex* __restrict__, int* __restrict__, int* __restrict__);
__device__ void getComplexInverse (const int, const int, cuDoubleComplex* __restrict__, const int* __restrict__,
										cuDoubleComplex* __restrict__);
__device__ void getComplexInverseHessenberg (const int, const int, cuDoubleComplex* __restrict__, int* __restrict__,
												cuDoubleComplex* __restrict__);

#endif