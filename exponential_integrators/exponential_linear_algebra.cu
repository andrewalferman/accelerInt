/* exponential_linear_algebra.c
 * Implementation of various linear algebra functions needed in the exponential integrators
 * \file exponential_linear_algebra.c
 *
 * \author Nicholas Curtis
 * \date 03/09/2015
 *
 */

#include "exponential_linear_algebra.cuh"

///////////////////////////////////////////////////////////////////////////////

/** Matrix-vector multiplication of a matrix sized MxM and a vector Mx1
 * 
 * \param[in]		m 		size of the matrix
 * \param[in]		A		matrix of size MxM
 * \param[in]		V		vector of size Mx1
 * \param[out]		Av		vector that is A * v
 */
__device__
void matvec_m_by_m (const int m, const double * const __restrict__ A,
						const double * const __restrict__ V,
						double * const __restrict__ Av) {
	//for each row
	for (int i = 0; i < m; ++i) {
		Av[INDEX(i)] = A[INDEX(i)] * V[INDEX(0)];
		
		//go across a row of A, multiplying by a column of phiHm
		for (int j = 1; j < m; ++j) {
			Av[INDEX(i)] += A[INDEX(j * STRIDE + i)] * V[INDEX(j)];
		}
	}
}

///////////////////////////////////////////////////////////////////////////////

/** Matrix-vector plus equals for a matrix of size MxM and vector of size Mx1
 * 
 *  That is, it returns (A + I) * v
 *
 * Performs inline matrix-vector multiplication (with unrolled loops) 
 * 
 * \param[in]		m 		size of the matrix
 * \param[in]		A		matrix of size MxM
 * \param[in]		V		vector of size Mx1
 * \param[out]		Av		vector that is (A + I) * v
 */
__device__ void matvec_m_by_m_plusequal (const int m, const double * const __restrict__ A, 
										 const double * const __restrict__ V, double * const __restrict__ Av)
{
	//for each row
	for (int i = 0; i < m; ++i) {
		Av[INDEX(i)] = A[INDEX(i)] * V[INDEX(0)];
		
		//go across a row of A, multiplying by a column of phiHm
		for (int j = 1; j < m; ++j) {
			Av[INDEX(i)] += A[INDEX(j * STRIDE + i)] * V[INDEX(j)];
		}

		Av[INDEX(i)] += V[INDEX(i)];
	}
}

/** Matrix-vector multiplication of a matrix sized NSPxM and a vector of size Mx1 scaled by a specified factor
 * 
 * Performs inline matrix-vector multiplication (with unrolled loops)
 * 
 * \param[in]		m 		size of the matrix
 * \param[in]		scale 	a number to scale the multplication by
 * \param[in]		A		matrix
 * \param[in]		V		the vector
 * \param[out]		Av		vector that is A * V
 */
__device__
void matvec_n_by_m_scale (const int m, const double scale,
						  const double * const __restrict__ A,
						  const double * const __restrict__ V,
						  double * const __restrict__ Av) {
	//for each row
	#pragma unroll
	for (int i = 0; i < NSP; ++i) {
		Av[INDEX(i)] = A[INDEX(i)] * V[INDEX(0)];
		
		//go across a row of A, multiplying by a column of phiHm
		for (int j = 1; j < m; ++j) {
			Av[INDEX(i)] += A[INDEX(j * NSP + i)] * V[INDEX(j)];
		}

		Av[INDEX(i)] *= scale;
	}
}


/** Matrix-vector multiplication of a matrix sized NSPxM and a vector of size Mx1 scaled by a specified factor
 *
 *  Computes the following:
 *  Av1 = A * V1 * scale[0]
 *  Av2 = A * V2 * scale[1]
 *  Av3 = A * V3 * scale[2] + V4 + V5
 * 
 * Performs inline matrix-vector multiplication (with unrolled loops)
 * 
 * \param[in]		m 		size of the matrix
 * \param[in]		scale 	a list of numbers to scale the multplication by
 * \param[in]		A		matrix
 * \param[in]		V		a list of 5 pointers corresponding to V1, V2, V3, V4, V5
 * \param[out]		Av		a list of 3 pointers corresponding to Av1, Av2, Av3
 */
__device__
void matvec_n_by_m_scale_special (const int m, const double * __restrict__ scale,
								  const double * __restrict__ A,
								  double * const __restrict__ * V,
								  double * __restrict__ * Av) {
	//for each row
	#pragma unroll
	for (int i = 0; i < NSP; ++i) {
		Av[0][INDEX(i)] = A[INDEX(i)] * V[0][INDEX(0)];
		Av[1][INDEX(i)] = A[INDEX(i)] * V[1][INDEX(0)];
		Av[2][INDEX(i)] = A[INDEX(i)] * V[2][INDEX(0)];
		
		//go across a row of A, multiplying by a column of phiHm
		for (int j = 1; j < m; ++j) {
			Av[0][INDEX(i)] += A[INDEX(j * NSP + i)] * V[0][INDEX(j)];
			Av[1][INDEX(i)] += A[INDEX(j * NSP + i)] * V[1][INDEX(j)];
			Av[2][INDEX(i)] += A[INDEX(j * NSP + i)] * V[2][INDEX(j)];
		}

		Av[0][INDEX(i)] *= scale[0];
		Av[1][INDEX(i)] *= scale[1];
		Av[2][INDEX(i)]  = scale[2] * Av[2][INDEX(i)] + V[3][INDEX(i)] + V[4][INDEX(i)];
	}
}

/** Matrix-vector multiplication of a matrix sized NSPxM and a vector of size Mx1 scaled by a specified factor
 *
 *  Computes the following:
 *  Av1 = A * V1 * scale[0]
 *  Av2 = A * V2 * scale[1]
 * 
 * Performs inline matrix-vector multiplication (with unrolled loops)
 * 
 * \param[in]		m 		size of the matrix
 * \param[in]		scale 	a list of numbers to scale the multplication by
 * \param[in]		A		matrix
 * \param[in]		V		a list of 2 pointers corresponding to V1, V2
 * \param[out]		Av		a list of 2 pointers corresponding to Av1, Av2
 */
__device__
void matvec_n_by_m_scale_special2 (const int m, const double* __restrict__ scale, const double* __restrict__ A,
										double* const __restrict__ * V, double* __restrict__ * Av) {
	//for each row
	#pragma unroll
	for (int i = 0; i < NSP; ++i) {
		Av[0][INDEX(i)] = A[INDEX(i)] * V[0][INDEX(0)];
		Av[1][INDEX(i)] = A[INDEX(i)] * V[1][INDEX(0)];
		
		//go across a row of A, multiplying by a column of phiHm
		for (int j = 1; j < m; ++j) {
			Av[0][INDEX(i)] += A[INDEX(j * NSP + i)] * V[0][INDEX(j)];
			Av[1][INDEX(i)] += A[INDEX(j * NSP + i)] * V[1][INDEX(j)];
		}

		Av[0][INDEX(i)] *= scale[0];
		Av[1][INDEX(i)] *= scale[1];
	}
}

///////////////////////////////////////////////////////////////////////////////

/** Matrix-vector multiplication of a matrix sized NSPxM and a vector of size Mx1 scaled by a specified factor and added to another vector
 * 
 * Performs inline matrix-vector multiplication (with unrolled loops)
 * 
 * \param[in]		m 		size of the matrix
 * \param[in]		scale 	a number to scale the multplication by
 * \param[in]		add 	the vector to add to the result
 * \param[in]		A		matrix
 * \param[in]		V		the vector
 * \param[out]		Av		vector that is A * V + add
 */
__device__
void matvec_n_by_m_scale_add (const int m, const double scale,
								const double* __restrict__ A, const double* __restrict__ V,
								double* __restrict__ Av, const double* __restrict__ add) {
	//for each row
	#pragma unroll
	for (int i = 0; i < NSP; ++i) {
		Av[INDEX(i)] = A[INDEX(i)] * V[INDEX(0)];
		
		//go across a row of A, multiplying by a column of phiHm
		for (int j = 1; j < m; ++j) {
			Av[INDEX(i)] += A[INDEX(j * NSP + i)] * V[INDEX(j)];
		}

		Av[INDEX(i)] = Av[INDEX(i)] * scale + add[INDEX(i)];
	}
}

///////////////////////////////////////////////////////////////////////////////

/** Matrix-vector multiplication of a matrix sized NSPxM and a vector of size Mx1 scaled by a specified factor and adds and subtracts the specified vectors
 *  note, the addition is twice the specified vector
 * 
 * Performs inline matrix-vector multiplication (with unrolled loops)
 * 
 * \param[in]		m 		size of the matrix
 * \param[in]		scale 	a number to scale the multplication by
 * \param[in]		add 	the vector to add to the result
 * \param[in]		sub 	the vector to subtract from the result
 * \param[in]		A		matrix
 * \param[in]		V		the vector
 * \param[out]		Av		vector that is scale * A * V + 2 * add - sub
 */
__device__
void matvec_n_by_m_scale_add_subtract (const int m, const double scale,
										const double* __restrict__ A, const double* V,
										double* __restrict__ Av, const double* __restrict__ add,
										const double* __restrict__ sub) {
	//for each row
	#pragma unroll
	for (int i = 0; i < NSP; ++i) {
		Av[INDEX(i)] = A[INDEX(i)] * V[INDEX(0)];
		
		//go across a row of A, multiplying by a column of phiHm
		#pragma unroll
		for (int j = 1; j < m; ++j) {
			Av[INDEX(i)] += A[INDEX(j * NSP + i)] * V[INDEX(j)];
		}

		Av[INDEX(i)] = Av[INDEX(i)] * scale + 2.0 * add[INDEX(i)] - sub[INDEX(i)];
	}
}

///////////////////////////////////////////////////////////////////////////////

/** Get scaling for weighted norm
 * 
 * \param[in]		y0		values at current timestep
 * \param[in]		y1		values at next timestep
 * \param[out]	sc	array of scaling values
 */
__device__
void scale (const double* __restrict__ y0, const double* __restrict__ y1, double* __restrict__ sc) {
	#pragma unroll
	for (int i = 0; i < NSP; ++i) {
		sc[INDEX(i)] = 1.0 / (ATOL + fmax(fabs(y0[INDEX(i)]), fabs(y1[INDEX(i)])) * RTOL);
	}
}

///////////////////////////////////////////////////////////////////////////////

/** Get scaling for weighted norm for the initial timestep (used in krylov process)
 * 
 * \param[in]		y0		values at current timestep
 * \param[out]	sc	array of scaling values
 */
__device__
void scale_init (const double* __restrict__ y0, double* __restrict__ sc) {
	#pragma unroll
	for (int i = 0; i < NSP; ++i) {
		sc[INDEX(i)] = 1.0 / (ATOL + fabs(y0[INDEX(i)]) * RTOL);
	}
}

///////////////////////////////////////////////////////////////////////////////

/** Perform weighted norm
 * 
 * \param[in]		nums	values to be normed
 * \param[in]		sc		scaling array for norm
 * \return			norm	weighted norm
 */
__device__
double sc_norm (const double* __restrict__ nums, const double* __restrict__ sc) {
	double norm = 0.0;
	
	#pragma unroll
	for (int i = 0; i < NSP; ++i) {
		norm += nums[INDEX(i)] * nums[INDEX(i)] * (sc[INDEX(i)] * sc[INDEX(i)]);
	}
	
	return sqrt(norm / ((double)NSP));
}

/** Computes and returns the two norm of a vector
 *
 *	\param[in]		v 		the vector
 */
__device__
double two_norm(const double* __restrict__ v)
{
	double norm = 0.0;
	#pragma unroll
	for (int i = 0; i < NSP; ++i) {
		norm += v[INDEX(i)] * v[INDEX(i)];
	}
	return sqrt(norm);
}

/** Normalize the input vector using a 2-norm
 * 
 * \param[in]		v		vector to be normalized
 * \param[out]		v_out	where to stick the normalized part of v (in a column)
 */
__device__
double normalize (const double* __restrict__ v, double* __restrict__ v_out) {
	
	double norm = two_norm(v);

	//unlikely to happen, if so, we still need to copy
	if (norm == 0.0)
		norm = 1.0;

	double m_norm = 1.0 / norm;

	#pragma unroll
	for (int i = 0; i < NSP; ++i) {
		v_out[INDEX(i)] = v[INDEX(i)] * m_norm;
	}
	return norm;
}


/** Performs the dot product of the w vector with the given Matrix
 * 
 * \param[in]		w   	the vector with with to dot
 * \param[in]		Vm		the subspace matrix
 * \out						the dot product of the specified vectors
 */
__device__
double dotproduct(const double* __restrict__ w, const double* __restrict__ Vm)
{
	double sum = 0;
	#pragma unroll
	for(int i = 0; i < NSP; i++)
	{
		sum += w[INDEX(i)] * Vm[INDEX(i)];
	}
	return sum;
}

/** Subtracts Vm scaled by s from w
 * 
 * \param[in]		s   	the scale multiplier to use
 * \param[in]		Vm		the subspace matrix
 * \param[out]		w 		the vector to subtract from
 */
__device__ void scale_subtract(const double s, const double* __restrict__ Vm, double* __restrict__ w)
{
	#pragma unroll
	for (int i = 0; i < NSP; i++)
	{
		w[INDEX(i)] -= s * Vm[INDEX(i)];
	}
}

/** Sets Vm to s * w
 * 
 * \param[in]		stride 	number of columns in Vm
 * \param[in]		s   	the scale multiplier to use
 * \param[in]		w 		the vector to use as a base
 * \param[out]		Vm		the subspace matrix to set
 */
__device__ void scale_mult(const double s, const double* __restrict__ w, double* __restrict__ Vm)
{
	#pragma unroll
	for (int i = 0; i < NSP; i++)
	{
		Vm[INDEX(i)] = w[INDEX(i)] * s;
	}
}