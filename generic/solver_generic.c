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

#include "header.h"
#include "solver.h"

#ifdef STIFF_METRICS
#include "stiffnessmetrics.h"
#include "math.h"
#endif

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

        #ifdef STIFF_METRICS
        double CEM = 0.0;
        double CSP;
        double stiffratio;
        double stiffindicator;
        int M;
        calculatemetrics(y_local, pr_local, &stiffratio, &stiffindicator, &CEM,
                        &CSP, &M, t, t_end);

        //StartTimer();
        double time0 = omp_get_wtime( );
        #endif

        // call integrator for one time step
        check_error(tid, integrate (t, t_end, pr_local, y_local));
        //double runtime = GetTimer();

        #ifdef STIFF_METRICS
        double runtime = omp_get_wtime( ) - time0;
        runtime /= 1000.0;

        int failflag = 0;
        #endif
        // update global array with integrated values
        for (int i = 0; i < NSP; i++)
        {
            #ifdef STIFF_METRICS
            #ifdef CHEM_UTILS_HEAD
            if (y_local[i] != y_local[i] || isinf(y_local[i]) || y_local[i] < (double) 0.0) {
              if (failflag == 0) {
                printf("Bad values:\n");
              }
              printf("%i,%.15e,",i,y_local[i]);
              //   printf("y_local:\n");
              //   for (int j = 0; j < NSP - 1; j++) {
              //     printf("%.15e,",y_local[j]);
              //   }
              //   printf("%.15e\n",y_local[NSP - 1]);
              //   printf("y_global:\n");
              //   for (int j = 0; j < NSP - 1; j++) {
              //     printf("%.15e,",y_global[tid + j * NUM]);
              //   }
              //   printf("%.15e\n",y_global[tid + (NSP - 1) * NUM]);
              }
              failflag = 1;
              printf("\n");
            }
            #else
            if (y_local[i] != y_local[i] || isinf(y_local[i])) {
              printf("%i,%.15e,",i,y_local[i]);
              //   printf("y_local:\n");
              //   for (int j = 0; j < NSP - 1; j++) {
              //     printf("%.15e,",y_local[j]);
              //   }
              //   printf("%.15e\n",y_local[NSP - 1]);
              //   printf("y_global:\n");
              //   for (int j = 0; j < NSP - 1; j++) {
              //     printf("%.15e,",y_global[tid + j * NUM]);
              //   }
              //   printf("%.15e\n",y_global[tid + (NSP - 1) * NUM]);
              }
              failflag = 1;
              printf("\n");
            }
            #endif
            // if (y_local[i] != y_local[i] || isinf(y_local[i]) || y_local[i] < (double) 0.0) {
            //   failflag = 1;
            // }
            #endif
            y_global[tid + i * NUM] = y_local[i];
        }
        // printf("%.15e\n", y_local[0]);
        #ifdef STIFF_METRICS
        if (failflag == 1) {
          runtime = -1;
        }

        // Print stiffness metrics and timing info
        printf("%i,%.15e,%.15e,%.15e,%.15e,%i,%.15e\n", tid, stiffratio, stiffindicator, CEM, CSP, M, runtime);
        #endif

    } //end tid loop
    //double runtime = GetTimer();
    //runtime /= 1000.0;
    //printf("Step: %.15e, Time: %.15e sec\n", t, runtime);
} // end intDriver

#ifdef GENERATE_DOCS
 }
#endif
