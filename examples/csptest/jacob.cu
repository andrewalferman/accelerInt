/**
 * \file
 * \brief A CUDA implementation of the van der Pol jacobian \f$\frac{\partial \dot{\vec{y}}}{\partial \vec{y}}\f$
 *
 * Implementes the CUDA evaluation of van der Pol Jacobian
 *
 */

#include "header.cuh"

#ifdef GENERATE_DOCS
//put this in the van der Pol namespace for documentation
namespace van_der_pol_cu {
#endif

/**
 * \brief An implementation of the van der Pol jacobian
 *
 * \param[in]           t               The current system time
 * \param[in]           mu              The van der Pol parameter
 * \param[in]           y               The state vector at time t
 * \param[out]          jac             The jacobian to populate
 * \param[in]           d_mem           The mechanism_memory struct.  In future versions, this will be used to access the \f$\mu\f$ parameter to have a consistent interface.
 *
 *  The Jacobian is in Column-major (Fortran) order.  As with dydt(), this function operates directly on
 *  global state vector and jacobian.  Hence we use the #INDEX macro defined in gpu_macros.cuh here.
 *  @see dydt()
 *
 */
__device__
void eval_jacob (const double t, const double eps, const double * __restrict__ y, double * __restrict__ jac, const mechanism_memory * __restrict__ d_mem)
{
    //Note, to reach index [i, j] of the Jacobian, we multiply `i` by NSP, the size of the first dimension of Jacobian and add j, i.e.:
    //jac[i, j] -> jac[i * NSP + j]
    jac[INDEX(0)]  = -1.0 / ( eps * eps * eps );
    jac[INDEX(4)]  = 1.0 / ( eps * eps * eps * ( 1.0 + y[INDEX(1)] ) * ( 1.0 + y[INDEX(1)] ) )
             + ( y[INDEX(1)] - 1.0 ) / ( (y[INDEX(1)] + 1.0) * (y[INDEX(1)] + 1.0) * (y[INDEX(1)] + 1.0) );
    jac[INDEX(8)]  = ( y[INDEX(2)] - 1.0 ) / ( (y[INDEX(2)] + 1.0) * (y[INDEX(2)] + 1.0) * (y[INDEX(2)] + 1.0) );
    jac[INDEX(12)] = ( y[INDEX(3)] - 1.0 ) / ( (y[INDEX(3)] + 1.0) * (y[INDEX(3)] + 1.0) * (y[INDEX(3)] + 1.0) );

    // dg1/dyi
  	/*
    jac[1]  = 0.0;
  	jac[5]  = -eps2;
  	jac[9]  = ( eps2 * y32 ) + ( ( y[2] - 1.0 ) * y33 );
  	jac[13] = ( y[3] - 1.0 ) * y43;
  	*/
    jac[INDEX(1)]  = 0.0;
    jac[INDEX(5)]  = -1.0 / ( eps * eps );
    jac[INDEX(9)]  = 1.0 / ( eps * eps * ( 1.0 + y[INDEX(2)] ) * ( 1.0 + y[INDEX(2)] ) )
             + ( y[INDEX(2)] - 1.0 ) / ( (y[INDEX(2)] + 1.0) * (y[INDEX(2)] + 1.0) * (y[INDEX(2)] + 1.0) );
    jac[INDEX(13)] = ( y[INDEX(3)] - 1.0 ) / ( (y[INDEX(3)] + 1.0) * (y[INDEX(3)] + 1.0) * (y[INDEX(3)] + 1.0) );

    // dg2/dyi
  	/*
    jac[2]  = 0.0;
  	jac[6]  = 0.0;
  	jac[10] = -eps1;
  	jac[14] = ( eps1 * y42 ) + ( ( y[3] - 1.0 ) * y43 );
  	*/
    jac[INDEX(2)]  = 0.0;
    jac[INDEX(6)]  = 0.0;
    jac[INDEX(10)] = -1.0 / eps;
    jac[INDEX(14)] = 1.0 / ( eps * ( 1.0 + y[INDEX(3)] ) * ( 1.0 + y[INDEX(3)] ) )
             + ( y[INDEX(3)] - 1.0 ) / ( (y[INDEX(3)] + 1.0) * (y[INDEX(3)] + 1.0) * (y[INDEX(3)] + 1.0) );

    // dg3/dyi
  	jac[INDEX(3)]  = 0.0;
  	jac[INDEX(7)]  = 0.0;
  	jac[INDEX(11)] = 0.0;
  	jac[INDEX(15)] = -1.0;
}

#ifdef GENERATE_DOCS
}
#endif
