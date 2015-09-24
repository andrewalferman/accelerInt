/** 
 * \file krylov.c
 *
 * \author Nicholas J. Curtis
 * \date 09/02/2014
 *
 * A krylov subspace integrator using the EXP4 method
 * based on the work of Niesen and Wright (2012)
 * 
 * NOTE: all matricies stored in column major format!
 * 
 */
 
/** Include common code. */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include <string.h>

#include "header.h"
#include "dydt.h"
#include "jacob.h"
#include "arnoldi.h"
#include "exp4_props.h"
#include "exponential_linear_algebra.h"
#include "solver_init.h"
#include "sparse_multiplier.h"

///////////////////////////////////////////////////////////////////////////////

/** 4th-order exponential integrator function w/ adaptive Kyrlov subspace approximation
 * 
 * 
 */
void integrate (const double t_start, const double t_end, const double pr, double* y) {
	
	//initial time
	double h = 1.0e-8;
	double h_max, h_new;

	double err_old = 1.0;
	double h_old = h;
	
	bool reject = false;

	double t = t_start;

	// get scaling for weighted norm
	double sc[NN];
	scale_init(y, sc);

#ifdef LOG_KRYLOV_AND_STEPSIZES
	//file for krylov logging
	FILE *logFile;
	//open and clear
	const char* f_name = solver_name();
	int len = strlen(f_name);
	char out_name[len + 17];
	sprintf(out_name, "log/%s-kry-log.txt", f_name);
	logFile = fopen(out_name, "w");

	char out_reject_name[len + 23];
	sprintf(out_reject_name, "log/%s-kry-reject.txt", f_name);	
	//file for krylov logging
	FILE *rFile;
	//open and clear
	rFile = fopen(out_reject_name, "a");
#endif
	
	double beta = 0;
	while ((t < t_end) && (t + h > t)) {
		//initial krylov subspace sizes
		int m, m1, m2;

		// temporary arrays
		double temp[NN];
		double f_temp[NN];
		double y1[NN];
		
		h_max = t_end - t;
	
		// source vector
		double fy[NN];
		dydt (t, pr, y, fy);

		// Jacobian matrix
		double A[NN * NN] = {0.0};
		eval_jacob (t, pr, y, A);

		double Hm[STRIDE * STRIDE] = {0.0};
		double Vm[NN * STRIDE] = {0.0};
		double phiHm[STRIDE * STRIDE] = {0.0};
		double err = 0.0;

		do
		{
			//do arnoldi
			if (arnoldi(&m, 1.0 / 3.0, P, h, A, fy, sc, &beta, Vm, Hm, phiHm) >= M_MAX)
			{
				h /= 3;
				continue;
			}

			// k1
			double k1[NN];
			//k1 is partially in the first column of phiHm
			//k1 = beta * Vm * phiHm(:, 1)
			matvec_n_by_m_scale(m, beta, Vm, phiHm, k1);
		
			// k2
			double k2[NN];
			//computing phi(2h * A)
			matvec_m_by_m (m, phiHm, phiHm, temp);
			//note: f_temp will contain hm * phi * phi * e1 for later use
			matvec_m_by_m (m, Hm, temp, f_temp);
			matvec_n_by_m_scale_add(m, beta * (h / 6.0), Vm, f_temp, k2, k1);
		
			// k3
			double k3[NN];
			//use the stored hm * phi * phi * e1 to get phi(3h * A)
			matvec_m_by_m (m, phiHm, f_temp, temp);
			matvec_m_by_m (m, Hm, temp, f_temp);
			matvec_n_by_m_scale_add_subtract(m, beta * (h * h / 27.0), Vm, f_temp, k3, k2, k1);
				
			// d4
			double k4[NN];
			#pragma unroll
			for (int i = 0; i < NN; ++i) {
				// f4
				f_temp[i] = h * ((-7.0 / 300.0) * k1[i] + (97.0 / 150.0) * k2[i] - (37.0 / 300.0) * k3[i]);
			
				k4[i] = y[i] + f_temp[i];
			}
			
			dydt (t, pr, k4, temp);
			sparse_multiplier (A, f_temp, k4);
		
			#pragma unroll
			for (int i = 0; i < NN; ++i) {
				k4[i] = temp[i] - fy[i] - k4[i];
			}

			//do arnoldi
			arnoldi(&m1, 1.0 / 3.0, P, h, A, k4, sc, &beta, Vm, Hm, phiHm);
			//k4 is partially in the m'th column of phiHm
			matvec_n_by_m_scale(m1, beta, Vm, phiHm, k4);
		
			// k5
			double k5[NN];
			//computing phi(2h * A)
			matvec_m_by_m (m1, phiHm, phiHm, temp);
			//note: f_temp will contain hm * phi * phi * e1 for later use
			matvec_m_by_m (m1, Hm, temp, f_temp);
			matvec_n_by_m_scale_add(m1, beta * (h / 6.0), Vm, f_temp, k5, k4);
				
			// k6
			double k6[NN];
			//use the stored hm * phi * phi * e1 to get phi(3h * A)
			matvec_m_by_m (m1, phiHm, f_temp, temp);
			matvec_m_by_m (m1, Hm, temp, f_temp);
			matvec_n_by_m_scale_add_subtract(m1, beta * (h * h / 27.0), Vm, f_temp, k6, k5, k4);
				
			// k7
			double k7[NN];
			#pragma unroll
			for (int i = 0; i < NN; ++i) {
				// f7
				f_temp[i] = h * ((59.0 / 300.0) * k1[i] - (7.0 / 75.0) * k2[i] + (269.0 / 300.0) * k3[i] + (2.0 / 3.0) * (k4[i] + k5[i] + k6[i]));
			
				k7[i] = y[i] + f_temp[i];
			}
		
			dydt (t, pr, k7, temp);
			sparse_multiplier (A, f_temp, k7);
		
			#pragma unroll
			for (int i = 0; i < NN; ++i) {
				k7[i] = temp[i] - fy[i] - k7[i];
			}
		
			arnoldi(&m2, 1.0 / 3.0, P, h, A, k7, sc, &beta, Vm, Hm, phiHm);
			//k7 is partially in the m'th column of phiHm
			matvec_n_by_m_scale(m2, beta / (h / 3.0), Vm, &phiHm[m2 * STRIDE], k7);
					
			// y_n+1
			#pragma unroll
			for (int i = 0; i < NN; ++i) {
				y1[i] = y[i] + h * (k3[i] + k4[i] - (4.0 / 3.0) * k5[i] + k6[i] + (1.0 / 6.0) * k7[i]);
			}
			
			scale (y, y1, f_temp);	
			
			///////////////////
			// calculate errors
			///////////////////
		
			// error of embedded order 3 method
			#pragma unroll
			for (int i = 0; i < NN; ++i) {
				temp[i] = k3[i] - (2.0 / 3.0) * k5[i] + 0.5 * (k6[i] + k7[i] - k4[i]) - (y1[i] - y[i]) / h;
			}	
			err = h * sc_norm(temp, f_temp);
			
			// error of embedded W method
			#pragma unroll
			for (int i = 0; i < NN; ++i) {
				temp[i] = -k1[i] + 2.0 * k2[i] - k4[i] + k7[i] - (y1[i] - y[i]) / h;
			}
			//double err_W = h * sc_norm(temp, sc);
			err = fmax(EPS, fmin(err, h * sc_norm(temp, f_temp)));
			
			// classical step size calculation
			h_new = pow(err, -1.0 / ORD);	
			
			if (err <= 1.0) {

				#ifdef LOG_KRYLOV_AND_STEPSIZES
					fprintf (logFile, "%.15le\t%.15le\t%.15le\t%d\t%d\t%d\n", t, h, err, m, m1, m2);
	  			#endif

				memcpy(sc, f_temp, NN * sizeof(double));
				
				// minimum of classical and Gustafsson step size prediction
				h_new = fmin(h_new, (h / h_old) * pow((err_old / (err * err)), (1.0 / ORD)));
				
				// limit to 0.2 <= (h_new/8) <= 8.0
				h_new = h * fmax(fmin(0.9 * h_new, 8.0), 0.2);
				h_new = fmin(h_new, h_max);
				
				// update y and t
				#pragma unroll
				for (int i = 0; i < NN; ++i) {
					y[i] = y1[i];
				}
				
				t += h;
				
				// store time step and error
				err_old = fmax(1.0e-2, err);
				h_old = h;
				
				// check if last step rejected
				if (reject) {
					reject = false;
					h_new = fmin(h, h_new);
				}
				h = fmin(h_new, fabs(t_end - t));
							
			} else {

				#ifdef LOG_KRYLOV_AND_STEPSIZES
					fprintf (rFile, "%.15le\t%.15le\t%.15le\t%d\t%d\t%d\n", t, h, err, m, m1, m2);
	  			#endif

				// limit to 0.2 <= (h_new/8) <= 8.0
				h_new = h * fmax(fmin(0.9 * h_new, 8.0), 0.2);
				h_new = fmin(h_new, h_max);

		
				reject = true;
				h = fmin(h, h_new);
			}
		} while(err >= 1.0);

	} // end while

	#ifdef LOG_KRYLOV_AND_STEPSIZES
		fclose(logFile);
		fclose(rFile);
	#endif
	
}