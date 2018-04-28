/**
 * \file
 * \brief An implementation of the van der Pol jacobian \f$\frac{\partial \dot{\vec{y}}}{\partial \vec{y}}\f$
 *
 * Implementes the evaluation of van der Pol Jacobian
 *
 */

#include "header.h"

#ifdef GENERATE_DOCS
//put this in the van der Pol namespace for documentation
namespace van_der_pol {
#endif

/**
 * \brief An implementation of the van der Pol jacobian
 *
 * \param[in]           t               The current system time
 * \param[in]           mu              The van der Pol parameter
 * \param[in]           y               The state vector at time t
 * \param[out]          jac             The jacobian to populate
 *
 *  The Jacobian is in a local Column-major (Fortran) order.  As with dydt(), this function operates on local
 *  copies of the global state vector and jacobian.  Hence simple linear indexing can be used here.
 *  @see solver_generic.c
 *
 */
void eval_jacob (const double t, const double eps, const double * __restrict__ y, double * __restrict__ jac)
{
    //Note, to reach index [i, j] of the Jacobian, we multiply `i` by NSP, the size of the first dimension of Jacobian and add j, i.e.:
    //jac[i, j] -> jac[i * NSP + j]
    jac[0]  = -1.0 / ( eps * eps * eps );
    jac[4]  = 1.0 / ( eps * eps * eps * ( 1.0 + y[1] ) * ( 1.0 + y[1] ) )
             + ( y[1] - 1.0 ) / ( (y[1] + 1.0) * (y[1] + 1.0) * (y[1] + 1.0) );
    jac[8]  = ( y[2] - 1.0 ) / ( (y[2] + 1.0) * (y[2] + 1.0) * (y[2] + 1.0) );
    jac[12] = ( y[3] - 1.0 ) / ( (y[3] + 1.0) * (y[3] + 1.0) * (y[3] + 1.0) );

    // dg1/dyi
  	/*
    jac[1]  = 0.0;
  	jac[5]  = -eps2;
  	jac[9]  = ( eps2 * y32 ) + ( ( y[2] - 1.0 ) * y33 );
  	jac[13] = ( y[3] - 1.0 ) * y43;
  	*/
    jac[1]  = 0.0;
    jac[5]  = -1.0 / ( eps * eps );
    jac[9]  = 1.0 / ( eps * eps * ( 1.0 + y[2] ) * ( 1.0 + y[2] ) )
             + ( y[2] - 1.0 ) / ( (y[2] + 1.0) * (y[2] + 1.0) * (y[2] + 1.0) );
    jac[13] = ( y[3] - 1.0 ) / ( (y[3] + 1.0) * (y[3] + 1.0) * (y[3] + 1.0) );

    // dg2/dyi
  	/*
    jac[2]  = 0.0;
  	jac[6]  = 0.0;
  	jac[10] = -eps1;
  	jac[14] = ( eps1 * y42 ) + ( ( y[3] - 1.0 ) * y43 );
  	*/
    jac[2]  = 0.0;
    jac[6]  = 0.0;
    jac[10] = -1.0 / eps;
    jac[14] = 1.0 / ( eps * ( 1.0 + y[3] ) * ( 1.0 + y[3] ) )
             + ( y[3] - 1.0 ) / ( (y[3] + 1.0) * (y[3] + 1.0) * (y[3] + 1.0) );

    // dg3/dyi
  	jac[3]  = 0.0;
  	jac[7]  = 0.0;
  	jac[11] = 0.0;
  	jac[15] = -1.0;
}

#ifdef GENERATE_DOCS
}
#endif
