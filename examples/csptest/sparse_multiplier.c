/**
 * \file
 * \brief Implementation for Jacobian vector multiplication, used in exponential integrators
 *
 */

#include "sparse_multiplier.h"

#ifdef GENERATE_DOCS
//put this in the van der Pol namespace for documentation
namespace van_der_pol {
#endif


/**
 * \brief Implements Jacobian \ vector multiplication in sparse (or unrolled) form
 * \param[in]           A           The (NSP x NSP) Jacobian matrix, see eval_jacob() for details on layout
 * \param[in]           Vm          The (NSP x 1) vector to multiply by
 * \param[out]          w           The (NSP x 1) vector to store the result in, \f$w := A * Vm\f$
 */
void sparse_multiplier(const double * A, const double * Vm, double* w) {
  w[0] =  A[0] * Vm[0] +  A[4] * Vm[1] + A[8] * Vm[2] + A[12] * Vm[3];
  w[1] =  A[1] * Vm[0] +  A[5] * Vm[1] + A[9] * Vm[2] + A[13] * Vm[3];
  w[2] =  A[2] * Vm[0] +  A[6] * Vm[1] + A[10] * Vm[2] + A[14] * Vm[3];
  w[3] =  A[3] * Vm[0] +  A[7] * Vm[1] + A[11] * Vm[2] + A[15] * Vm[3];
}


#ifdef GENERATE_DOCS
}
#endif
