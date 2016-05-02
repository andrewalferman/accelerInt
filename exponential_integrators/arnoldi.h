/* krylov.h
 * Implementation of the arnoldi iteration methods
 * \file arnoldi.h
 *
 * \author Nicholas Curtis
 * \date 03/09/2015
 *
 */

#ifndef ARNOLDI_H
#define ARNOLDI_H

#include <string.h>

#include "header.h"
#include "phiAHessenberg.h"
#include "exponential_linear_algebra.h"
#include "sparse_multiplier.h"
#include "solver_options.h"
#include "solver_props.h"

static int index_list[23] = {1, 2, 3, 4, 5, 6, 7, 9, 11, 14, 17, 21, 27, 34, 42, 53, 67, 84, 106, 133, 167, 211, 265};

///////////////////////////////////////////////////////////////////////////////

/** Matrix-vector multiplication of a matrix sized MxM and a vector Mx1
 * 
 * Performs inline matrix-vector multiplication (with unrolled loops) 
 * 
 * \param[in, out]		m 		in - the starting size of the matrix, out - the ending size of the matrix
 * \param[in]			scale	the value to scale the timestep by
 * \param[in]			p		the order of the maximum phi function needed
 * \param[in]			h		the timestep
 * \param[in]			A 		the jacobian matrix
 * \param[in]  			v 		the vector to use for the krylov subspace
 * \param[in] 			sc 		the error scaling vector
 * \param[out] 			beta 	the norm of the v vector
 * \param[out]			Vm 		the arnoldi basis matrix
 * \param[out]			Hm 		the constructed Hm matrix, used in actual expoentials
 * \param[out] 			phiHm   the exponential matrix computed from h * scale * Hm
 */
static inline
int arnoldi(const double scale, const int p, const double h, const double* A, const double* v, const double* sc, double* beta, double* Vm, double* Hm, double* phiHm)
{
	//the temporary work array
	double w[NSP];

	//first place A*fy in the Vm matrix
	*beta = normalize(v, Vm);

	double store = 0;
	int index = 0;
	int j = 0;
	double err = 2.0;

	while(err > 1.0)
	{
		
		for (; j < index_list[index] && j + p < STRIDE; j++)
		{
			sparse_multiplier(A, &Vm[j * NSP], w);
			for (int i = 0; i <= j; i++)
			{
				Hm[j * STRIDE + i] = dotproduct(w, &Vm[i * NSP]);
				scale_subtract(Hm[j * STRIDE + i], &Vm[i * NSP], w);
			}
			Hm[j * STRIDE + j + 1] = two_norm(w);
			if (fabs(Hm[j * STRIDE + j + 1]) < ATOL)
			{
				//happy breakdown
				j++;
				break;
			}
			scale_mult(1.0 / Hm[j * STRIDE + j + 1], w, &Vm[(j + 1) * NSP]);
		}
		if (j + p >= STRIDE)
			return j;
		index++;
		//resize Hm to be mxm, and store Hm(m, m + 1) for later
		store = Hm[(j - 1) * STRIDE + j];
		Hm[(j - 1) * STRIDE + j] = 0.0;

		//0. fill potentially non-empty memory first
		memset(&Hm[j * STRIDE], 0, (j + 2) * sizeof(double));

		//get error
		//1. Construct augmented Hm (fill in identity matrix)
		Hm[(j) * STRIDE] = 1.0;
		
		for (int i = 1; i < p; i++)
		{
			//0. fill potentially non-empty memory first
			memset(&Hm[(j + i) * STRIDE], 0, (j + i + 2) * sizeof(double));
			Hm[(j + i) * STRIDE + (j + i - 1)] = 1.0;
		}

#ifdef RB43
		//2. Get phiHm
		int info = expAc_variable (j + p, Hm, h * scale, phiHm);
#elif EXP4
		//2. Get phiHm
		int info = phiAc_variable (j + p, Hm, h * scale, phiHm);
#endif
		if (info != 0)
		{
			return -info;
		}

		//3. Get error
		err = h * (*beta) * fabs(store * phiHm[(j) * STRIDE + (j) - 1]) * sc_norm(&Vm[(j) * NSP], sc);

		//restore Hm(m, m + 1)
		Hm[(j - 1) * STRIDE + j] = store;
		//restore real Hm (NOTE: untested in RB43)
		Hm[(j) * STRIDE] = 0.0;

	}

	return j;
}

#endif