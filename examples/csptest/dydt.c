/**
 * \file
 * \brief An implementation of the van der Pol right hand side (y' = f(y)) function.
 *
 * Implements the evaluation of the derivative of the state vector with respect to time, \f$\dot{\vec{y}}\f$
 *
 */

#include "header.h"

#ifdef GENERATE_DOCS
//put this in the van der Pol namespace for documentation
namespace van_der_pol {
#endif

/**
 * \brief An implementation of the RHS of the van der Pol equation
 * \param[in]        t         The current system time
 * \param[in]        mu        The van der Pol parameter
 * \param[in]        y         The state vector
 * \param[out]       dy        The output RHS (dydt) vector
 *
 * The `y` and `dy` vectors supplied here are local versions of the global state vectors.
 * They have been transformed from the global Column major (Fortan) ordering to a local 1-D vector
 * Hence the vector accesses can be done in a simple manner below, i.e. y[0] -> \f$y_1\f$, y[1] -> \f$y_2\f$, etc.
 * @see solver_generic.c
 */
void dydt (const double t, const double eps, const double * __restrict__ y, double * __restrict__ dy) {

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
  dy[0] = ( ( ( 1.0 / (eps * eps * eps) ) * ( -y[0] + ( y[1] / (1.0 + y[1]) ) )
          - ( y[1] / ( (y[1] + 1.0) * (y[1] + 1.0) ) ) )
          - ( y[2] / ( (y[2] + 1.0) * (y[2] + 1.0) ) ) )
          - ( y[3] / ( (y[3] + 1.0) * (y[3] + 1.0) ) );

  dy[1] = ( ( 1.0 / (eps * eps) ) * ( -y[1] + ( y[2] / (1.0 + y[2]) ) )
          - ( y[2] / ( (y[2] + 1.0) * (y[2] + 1.0) ) ) )
          - ( y[3] / ( (y[3] + 1.0) * (y[3] + 1.0) ) );

  dy[2] = ( 1.0 / eps ) * ( -y[2] + ( y[3] / (1.0 + y[3]) ) )
          - ( y[3] / ( (y[3] + 1.0) * (y[3] + 1.0) ) );

  dy[3] = -y[3];

} // end dydt


#ifdef GENERATE_DOCS
}
#endif
