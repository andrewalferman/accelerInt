/**
 * \file
 * \brief the generic integration driver for the CPU solvers
 *
 * \author Nicholas Curtis
 * \date 03/10/2015
 *
 * Modifications from Andrew Alferman
 * 11/14/2017
 *
 */
//#include <omp.h>

#include "header.h"
#include "solver.h"
//#include "timer.h"
#include "jacob.h"

/* DGEEV prototype */
extern void dgeev( char* jobvl, char* jobvr, int* n, double* a,
                int* lda, double* wr, double* wi, double* vl, int* ldvl,
                double* vr, int* ldvr, double* work, int* lwork, int* info );

/* Parameters */
#define N NSP
#define LDA N
#define LDVL N
#define LDVR N

#ifdef GENERATE_DOCS
 namespace generic {
#endif

/**
 * \brief Integration driver for the CPU integrators
 * \param[in]       NUM             The (non-padded) number of IVPs to integrate
 * \param[in]       t               The current system time
 * \param[in]       t_end           The IVP integration end time
 * \param[in]       pr_global       The system constant variable (pressures / densities)
 * \param[in,out]   y_global        The system state vectors at time t.
                                    Returns system state vectors at time t_end
 *
 * This is generic driver for CPU integrators
 */
void intDriver (const int NUM, const double t, const double t_end,
                const double *pr_global, double *y_global)
{
    //StartTimer();
    int tid;
    #pragma omp parallel for shared(y_global, pr_global) private(tid)
    for (tid = 0; tid < NUM; ++tid) {

        // local array with initial values
        double y_local[NSP];
        double pr_local = pr_global[tid];

        // load local array with initial values from global array

        for (int i = 0; i < NSP; i++)
        {
            y_local[i] = y_global[tid + i * NUM];
        }

        //StartTimer();
        double time0 = omp_get_wtime( );
        // call integrator for one time step
        check_error(tid, integrate (t, t_end, pr_local, y_local));
        //double runtime = GetTimer();
        double runtime = omp_get_wtime( ) - time0;
        runtime /= 1000.0;

        int failflag = 0;
        // update global array with integrated values
        for (int i = 0; i < NSP; i++)
        {
            y_global[tid + i * NUM] = y_local[i];
            if (isnan(y_local[i])) {
              failflag = 1;
            }
        }

        // if (tid == 296) {
        //   printf("Y Vector:\n");
        //   for (int i = 0; i < NSP; i++) {
        //     printf("%.15e\n", y_local[i]);
        //   }
        // }
        //

        if (failflag == 0) {
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
          //
          // if (tid == 296) {
          //   printf("Jacobian:\n");
          //   for (int i = 0; i < NSP*NSP; i++) {
          //     printf("%.15e\n", jac[i]);
          //   }
          //   printf("Hermitian:\n");
          //   for (int i = 0; i < NSP*NSP; i++) {
          //     printf("%.15e\n", hermitian[i]);
          //   }
          // }
          //
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

          // Print stiffness metrics and timing info
          printf("%i,%.15e,%.15e,%.15e,%.15e\n", tid, stiffratio, stiffindicator, CEM, runtime);
          // /* Print eigenvalues */
          // print_eigenvalues( "Eigenvalues", n, wr, wi );

          //Test print statement
          //printf("%.15e, %.15e\n", jacobian[0][0], jacobian[0][10])
        }
        else {
          printf("%i,-1.0,-1.0,-1.0,-1.0");
        }

    } //end tid loop
    //double runtime = GetTimer();
    //runtime /= 1000.0;
    //printf("Step: %.15e, Time: %.15e sec\n", t, runtime);
} // end intDriver

#ifdef GENERATE_DOCS
 }
#endif
