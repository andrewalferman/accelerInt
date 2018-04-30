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
#include "stiffnessmetrics.h"
#endif

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
        double CEM = 0.0;
        double CSP;
        double stiffratio;
        double stiffindicator;
        int M;
        calculatemetrics(y_local, pr_local, &stiffratio, &stiffindicator, &CEM,
                        &CSP, &M, t, t_end);
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
        // runtime /= 1000.0;
        //printf("Temp: %.15e, Time: %.15e sec\n", y_local[0], runtime);

        // update global array with integrated values and print output
        //printf("%i,", tid);

        int failflag = 0;
        #endif

        for (int i = 0; i < NSP; i++)
        {
            #ifdef STIFF_METRICS
            #ifdef CHEM_UTILS_HEAD
            if (y_local[i] != y_local[i] || isinf(y_local[i]) || y_local[i] < (double) 0.0) {
              if (failflag == 0) {
                printf("Bad values:\n");
              }
              printf("%i,%.15e,",i,y_local[i]);
              failflag = 1;
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

            #else
            if (y_local[i] != y_local[i] || isinf(y_local[i])) {
              if (failflag == 0) {
                printf("Bad values:\n");
              }
              printf("%i,%.15e,",i,y_local[i]);
              failflag = 1;
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
          printf("\n");
          runtime = -1;
        }

        // Print stiffness metrics and timing info
        #ifdef CHEM_UTILS_HEAD
        //printf("%i,%.15e,%.15e,%.15e,%.15e,%i,%.15e,%.15e\n", tid, stiffratio, stiffindicator, CEM, CSP, M, runtime,y_local[0]);
        printf("%i,%.15e,%.15e\n", tid, runtime, y_local[0]);
        #else
        //printf("%i,%.15e,%.15e,%.15e,%.15e,%i,%.15e,", tid, stiffratio, stiffindicator, CEM, CSP, M, runtime);
        printf("%i,%.15e,", tid, runtime);
        for (int i = 0; i < NSP; ++i) {
          printf("%.15e,", y_local[i]);
        }
        printf("\n");
        #endif
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
