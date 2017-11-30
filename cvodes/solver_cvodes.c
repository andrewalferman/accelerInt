/**
 * \file
 * \brief The integration driver for the CVODE solver
 *
 * \author Nicholas Curtis
 * \date 03/10/2015
 *
 * Modifications from Andrew Alferman
 * 11/14/2017
 *
 */

#include <omp.h>
#include <stdlib.h>

#include "header.h"
#include "solver.h"
#ifndef STIFF_METRICS
#include "timer.h"
#endif

/* CVODES INCLUDES */
#include "sundials/sundials_types.h"
#include "sundials/sundials_math.h"
#include "sundials/sundials_nvector.h"
#include "nvector/nvector_serial.h"
#include "cvodes/cvodes.h"
#include "cvodes/cvodes_lapack.h"
#include "cvodes_jac.h"
//#include "cvodes_dydt.h"
#include "jacob.h"

extern N_Vector *y_locals;
extern double* y_local_vectors;
extern void** integrators;

#ifdef STIFF_METRICS
/* DGEEV prototype */
extern void dgeev( char* jobvl, char* jobvr, int* n, double* a,
                int* lda, double* wr, double* wi, double* vl, int* ldvl,
                double* vr, int* ldvr, double* work, int* lwork, int* info );
#endif



/* Parameters */
#define N NSP
#define LDA N
#define LDVL N
#define LDVR N

#ifdef GENERATE_DOCS
namespace cvode {
#endif

/**
 * \brief Integration driver for the CPU integrators
 * \param[in]       NUM         the number of IVPs to solve
 * \param[in]       t           the current IVP time
 * \param[in]       t_end       the time to integrate the IVP to
 * \param[in]       pr_global   the pressure value for the IVPs
 * \param[in, out]  y_global    the state vectors
 *
 * The integration driver for the CVODEs solver
 */
void intDriver (const int NUM, const double t, const double t_end,
                const double *pr_global, double *y_global)
{
    int tid;
    double t_next;
    #pragma omp parallel for shared(y_global, pr_global, integrators, y_locals) private(tid, t_next)
    for (tid = 0; tid < NUM; ++tid) {
        int index = omp_get_thread_num();

        // local array with initial values
        N_Vector fill = y_locals[index];
        double pr_local = pr_global[tid];

        // load local array with initial values from global array
        double* y_local = NV_DATA_S(fill);

        for (int i = 0; i < NSP; i++)
        {
            y_local[i] = y_global[tid + i * NUM];
        }

        #ifdef STIFF_METRICS
        // Calculate the stiffness metrics
      	double jac[NSP*NSP];
      	eval_jacob(t_end, pr_global[tid], y_local, jac);
        // Get the Hermitian
        double hermitian[NSP*NSP];
        for (int i = 0; i < NSP; i++) {
          for (int j = 0; j < NSP; j++) {
            //jacobian[i][j] = jac[i * NSP + j];
            hermitian[i * NSP + j] = 0.5 * (jac[i * NSP + j] + jac[j * NSP + i]);
          }
        }

        // Get the eigenvalues of both matrices
        // Made 2 sets of variables in case dgeev messes them up
        int n = N, lda = LDA, ldvl = LDVL, ldvr = LDVR, info, lwork;
        int nh = N, ldah = LDA, ldvlh = LDVL, ldvrh = LDVR, infoh, lworkh;
        double wkopt;
        double* work;
        double wkopth;
        double* workh;
        /* Local arrays */
        double wr[N], wi[N], vl[LDVL*N], vr[LDVR*N];
        double xr[N], xi[N], ul[LDVL*N], ur[LDVR*N];

        /* Query and allocate the optimal workspace */
        // First, the Jacobian
        lwork = -1;
        dgeev( "Vectors", "Vectors", &n, jac, &lda, wr, wi, vl, &ldvl, vr, &ldvr,
         &wkopt, &lwork, &info );
        lwork = (int)wkopt;
        work = (double*)malloc( lwork*sizeof(double) );
        /* Solve eigenproblem */
        dgeev( "Vectors", "Vectors", &n, jac, &lda, wr, wi, vl, &ldvl, vr, &ldvr,
         work, &lwork, &info );
        /* Check for convergence */
        if( info > 0 ) {
                printf( "The algorithm failed to compute eigenvalues.\n" );
                exit( 1 );
        }
        // Next, the Hermitian
        lworkh = -1;
        dgeev( "Vectors", "Vectors", &nh, hermitian, &ldah, xr, xi, ul, &ldvlh, ur, &ldvrh,
         &wkopth, &lworkh, &infoh );
        lworkh = (int)wkopth;
        workh = (double*)malloc( lworkh*sizeof(double) );
        /* Solve eigenproblem */
        dgeev( "Vectors", "Vectors", &nh, hermitian, &ldah, xr, xi, ul, &ldvlh, ur, &ldvrh,
         workh, &lworkh, &infoh );
        /* Check for convergence */
        if( info > 0 ) {
                printf( "The algorithm failed to compute eigenvalues.\n" );
                exit( 1 );
        }

        // Get the minimum and maximum values of the eigenvalues
        double minjaceig = 1.0e10;
        double CEM = 0.0;
        double maxjaceig = 0.0;
        double minhereig = 1.0e10;
        double maxhereig = 0.0;
        for (int i = 0; i < NSP; i++) {
          if (abs(wr[i]) < minjaceig && abs(wr[i]) > 2.22045e-16) {
            minjaceig = abs(wr[i]);
          }
          if (wr[i] > CEM) {
            CEM = wr[i];
          }
          if (abs(wr[i]) > maxjaceig) {
            maxjaceig = abs(wr[i]);
          }
          if (xr[i] < minhereig) {
            minhereig = xr[i];
          }
          if (xr[i] > maxhereig) {
            maxhereig = xr[i];
          }
        }

        double stiffratio = maxjaceig / minjaceig;
        double stiffindicator = 0.5 * (minhereig + maxhereig);
        #endif

        //reinit this integrator for time t, w/ updated state
        int flag = CVodeReInit(integrators[index], t, fill);
        if (flag != CV_SUCCESS)
        {
            printf("Error reinitializing integrator for thread %d, code: %d\n", tid, flag);
            exit(flag);
        }

        //set user data to Pr
        flag = CVodeSetUserData(integrators[index], &pr_local);
        if (flag != CV_SUCCESS)
        {
            printf("Error setting user data for thread %d, code: %d\n", tid, flag);
            exit(flag);
        }

        //set end time
        flag = CVodeSetStopTime(integrators[index], t_end);
        if (flag != CV_SUCCESS)
        {
            printf("Error setting end time for thread %d, code: %d\n", tid, flag);
            exit(flag);
        }

        #ifdef STIFF_METRICS
        // Need to replace this with a threadsafe non-OMP method
        double time0 = omp_get_wtime( );
        #endif

        // call integrator for one time step
        flag = CVode(integrators[index], t_end, fill, &t_next, CV_NORMAL);
        if ((flag != CV_SUCCESS && flag != CV_TSTOP_RETURN) || t_next != t_end)
        {
            printf("Error on integration step for thread %d, code %d\n", tid, flag);
            exit(flag);
        }
        #ifdef STIFF_METRICS
        double runtime = omp_get_wtime( ) - time0;
        runtime /= 1000.0;
        //printf("Temp: %.15e, Time: %.15e sec\n", y_local[0], runtime);

        // update global array with integrated values and print output
        //printf("%i,", tid);

        int failflag = 0;
        #endif

        for (int i = 0; i < NSP; i++)
        {
            y_global[tid + i * NUM] = y_local[i];
            #ifdef STIFF_METRICS
            if (y_local[i] != y_local[i] || isinf(y_local[i]) || y_local[i] < (double) 0.0) {
              failflag = 1;
            }
            #endif
        }

        #ifdef STIFF_METRICS
        if (failflag == 1) {
          runtime = -1;
        }

        // Print stiffness metrics and timing info
        printf("%i,%.15e,%.15e,%.15e,%.15e\n", tid, stiffratio, stiffindicator, CEM, runtime);
        // /* Print eigenvalues */
        // print_eigenvalues( "Eigenvalues", n, wr, wi );

        //Test print statement
        //printf("%.15e, %.15e\n", jacobian[0][0], jacobian[0][10])
        #endif

    } // end tid loop

} // end intDriver

#ifdef GENERATE_DOCS
}
#endif
