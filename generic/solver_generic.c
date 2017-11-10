/**
 * \file
 * \brief the generic integration driver for the CPU solvers
 *
 * \author Nicholas Curtis
 * \date 03/10/2015
 *
 */
#include <omp.h>

#include "header.h"
#include "solver.h"
//#include "timer.h"

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
        //printf("%.15e,%.15e\n", y_global[tid], runtime);
        printf("%i,", tid);
        // update global array with integrated values

        for (int i = 0; i < NSP; i++)
        {
            y_global[tid + i * NUM] = y_local[i];
            printf("%.15e,", y_local[i]);
        }
        printf("%.15e\n", runtime);

    } //end tid loop
    //double runtime = GetTimer();
    //runtime /= 1000.0;
    //printf("Step: %.15e, Time: %.15e sec\n", t, runtime);
} // end intDriver

#ifdef GENERATE_DOCS
 }
#endif
