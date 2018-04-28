/**
 * \file
 * \brief A CUDA implementation of the van der Pol right hand side (y' = f(y)) function.
 *
 * Implements the CUDA evaluation of the derivative of the state vector with respect to time, \f$\dot{\vec{y}}\f$
 *
 */

#include "header.cuh"
#include "gpu_macros.cuh"

#ifdef GENERATE_DOCS
//put this in the van der Pol namespace for documentation
namespace van_der_pol_cu {
#endif

/**
 * \brief An implementation of the RHS of the van der Pol equation
 * \param[in]        t         The current system time
 * \param[in]        mu        The van der Pol parameter
 * \param[in]        y         The state vector
 * \param[out]       dy        The output RHS (dydt) vector
 * \param[in]        d_mem     The mechanism_memory struct.  In future versions, this will be used to access the \f$\mu\f$ parameter to have a consistent interface.
 *
 * The `y` and `dy` vectors supplied here are the global state vectors.
 * Hence the vector accesses must be done with the global thread ID.
 * The gpu_macros.cuh file defines some useful macros to simplify indexing
 * @see gpu_macros.cuh
 */
 __device__
void dydt (const double t, const double eps, const double * __restrict__ y, double * __restrict__ dy,
           const mechanism_memory * __restrict__ d_mem) {

     // Real eps1 = ONE / eps;
     // Real eps2 = ONE / ( eps * eps );
     //
     // Real y2 = ONE + y[1];
     //
     // Real y3 = ONE + y[2];
     // Real y32 = y3 * y3;
     //
     // Real y4 = ONE + y[3];
     // Real y42 = y4 * y4;

     /*
     ydot[0] = ( eps1 * eps2 ) * ( -y[0] + ( y[1] / y2 ) )
             - ( y[1] / ( y2 * y2 ) ) - ( y[2] / y32 )
             - ( y[3] / y42 );

     ydot[1] = eps2 * ( -y[1] + ( y[2] / y3 ) )
             - ( y[2] / y32 ) - ( y[3] / y42 );

     ydot[2] = eps1 * ( -y[2] + ( y[3] / y4 ) )
             - ( y[3] / y42 );
     */
     dy[INDEX(0)] = ( ( ( 1.0 / (eps * eps * eps) ) * ( -y[INDEX(0)] + ( y[INDEX(1)] / (1.0 + y[INDEX(1)]) ) )
             - ( y[INDEX(1)] / ( (y[INDEX(1)] + 1.0) * (y[INDEX(1)] + 1.0) ) ) )
             - ( y[INDEX(2)] / ( (y[INDEX(2)] + 1.0) * (y[INDEX(2)] + 1.0) ) ) )
             - ( y[INDEX(3)] / ( (y[INDEX(3)] + 1.0) * (y[INDEX(3)] + 1.0) ) );

     dy[INDEX(1)] = ( ( 1.0 / (eps * eps) ) * ( -y[INDEX(1)] + ( y[INDEX(2)] / (1.0 + y[INDEX(2)]) ) )
             - ( y[INDEX(2)] / ( (y[INDEX(2)] + 1.0) * (y[INDEX(2)] + 1.0) ) ) )
             - ( y[INDEX(3)] / ( (y[INDEX(3)] + 1.0) * (y[INDEX(3)] + 1.0) ) );

     dy[INDEX(2)] = ( 1.0 / eps ) * ( -y[INDEX(2)] + ( y[INDEX(3)] / (1.0 + y[INDEX(3)]) ) )
             - ( y[INDEX(3)] / ( (y[INDEX(3)] + 1.0) * (y[INDEX(3)] + 1.0) ) );

     dy[INDEX(3)] = -y[INDEX(3)];

} // end dydt


#ifdef GENERATE_DOCS
}
#endif
